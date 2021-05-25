#include <string>
#include <iostream>
#include <set>
#include <map>
#include <vector>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

using namespace amrex;

using std::vector;
using std::map;
using std::string;
using std::set;
using std::cerr;
using std::endl;
using std::cout;
using std::pair;
using std::ofstream;
using std::ostream;

struct MLloc
{
    MLloc() :
        amr_lev(-1), box_idx(-1), pt_idx(-1) {}
    MLloc(int lev,int box,int pt) :
        amr_lev(lev), box_idx(box), pt_idx(pt) {}
    int amr_lev, box_idx, pt_idx;
};

#ifdef BL_WRITE_BINARY
#include "TECIO.h"
#define SIZET INTEGER4
void
write_binary_tec_file(const std::string&          outfile,
                      const std::vector<string>&  names,
                      const Vector<int>&           faceData,
                      const Vector<MultiFab*>&     nodes,
                      const std::vector<MLloc>&   nodeMap,
                      const Vector<Vector<Real> >&  eltData,
                      const std::string&          title,
                      const Vector<int>&           idX);
#else
#define SIZET int
void
write_ascii_tec_file(const std::string&          outfile,
                     const std::vector<string>&  names,
                     const Vector<int>&           faceData,
                     const Vector<MultiFab*>&     nodes,
                     const std::vector<MLloc>&   nodeMap,
                     const Vector<Vector<Real> >&  eltData,
                     const std::string&          title,
                     const Vector<int>&           idX);
#endif

void
write_binary_mef_file(const std::string&          outfile,
                      const std::vector<string>&  names,
                      const Vector<int>&           faceData,
                      const Vector<MultiFab*>&     nodes,
                      const std::vector<MLloc>&   nodeMap,
                      const Vector<Vector<Real> >&  eltData,
                      const std::string&          title,
                      const Vector<int>&           idX);

int
ml_Nlevels(const std::string& infile);

void
read_ml_streamline_names(const std::string& infile,
                         Vector<int>&        comps,
                         vector<string>&    names);

void
read_ml_streamline_data(const std::string& infile,
                        Vector<int>&        strComps,
                        Vector<MultiFab*>&  data,
                        vector<string>&    names,
                        Vector<int>&        faceData,
                        int&               nElts,
                        Vector<Vector<Vector<int> > >& inside_nodes);
std::vector<MLloc>
build_nodeMap(const Vector<Vector<Vector<int> > >& inside_nodes);

Real
wedge_volume_int(const Vector<const MLloc*>& p,
                 int                        ptOnStr,
                 int                        comp,
                 const Vector<MultiFab*>&    data,
                 const Vector<int>& idX);

Real
wedge_surf_area(Vector<const MLloc*> p,
                const Vector<MultiFab*>& data,
                const Vector<int>& idX, int ptOnStr);

int get_nPts(const Vector<MultiFab*>& data);
int get_jlo(const Vector<MultiFab*>& data);

void
read_iso(const std::string& infile,
         FArrayBox&         nodes,
         Vector<int>&        faceData,
         int&               nElts,
         vector<string>&    names,
         string&            label);

Real
max_grad(const MLloc&            p,
         int                     comp,
         const Vector<MultiFab*>& data,
         const Vector<int>&       idX);

void
peak_val(const MLloc&            p,
         int                     pComp,
         Vector<int>&             sampleComps,
         Vector<Real>&            peakSamples,
         bool*                   peakOK,
         const Vector<MultiFab*>& data);

struct Edge
{
    Edge(int i1, int i2)
        {
            BL_ASSERT(i1!=i2);
            p = ( i1<i2  ?  pair<int,int>(i1,i2)  :  pair<int,int>(i2,i1)  );
        }
    Edge(const Edge& e)
        : p(e.p) {}
    pair<int,int> p;
};

struct EdgeLT
{
    bool operator()(const Edge& lhs, const Edge& rhs) const
        {
            if (lhs.p.first > rhs.p.first)
                return false;
            if (lhs.p.first == rhs.p.first)
                return lhs.p.second<rhs.p.second;
            else
                return true;
        }
};

typedef map<Edge,pair<int,int>,EdgeLT > NeighborMap;

ostream& operator<<(ostream& os, const NeighborMap& m)
{
    os << "NeighborMap: \n";
    for (NeighborMap::const_iterator it=m.begin(); it!=m.end(); ++it)
    {
        const Edge& e = it->first;
        const pair<int,int>& p=it->second;
        os << "  Edge(" << e.p.first << ", " << e.p.second << ") between elts: "
           << p.first << ", " << p.second << '\n';
    }
    return os;
}


static 
Vector<Vector<int> > 
buildNodeNeighbors(const Vector<int>& faceData,
                   int               nNodes)
{
    Vector<Vector<int> > nodeNeighbors(nNodes);
    int nElts = faceData.size()/3;
    for (int i=0; i<nElts; ++i)
    {
        int offset=i*3;
        for (int j=0; j<3; ++j)
        {
            Vector<int>& neighbors = nodeNeighbors[ faceData[offset+j]-1 ];
            neighbors.resize(neighbors.size()+1);
            neighbors[neighbors.size()-1] = i; // b/c of how this is done, we know the list will be unique
        }
    }

    Vector<Vector<int> > cellNeighbors(nElts);
    for (int i=0; i<nElts; ++i)
    {
        int offset=i*3;
        set<int> n;
        for (int j=0; j<3; ++j)
        {
            const Vector<int>& neighbors = nodeNeighbors[ faceData[offset+j]-1 ];
            for (int k=0; k<neighbors.size(); ++k)
            {
                int nc = neighbors[k];
                if (nc != i)
                    n.insert(nc);
            }
        }
        cellNeighbors[i].resize(n.size());
        int cnt=0;
        for (std::set<int>::const_iterator it=n.begin(); it!=n.end(); ++it)
            cellNeighbors[i][cnt++] = *it;
    }
    return cellNeighbors;
}


void
smoothVals(Vector<Real>&              newVals,
           const Vector<Real>&        vals,
           const Vector<Real>&        area,
           const Vector<Vector<int> >& nodeNeighbors)
{
    int nElts = vals.size();
    if (nElts!=newVals.size() || nElts!=area.size())
        Abort("Data wrong size");

    for (int i=0; i<nElts; ++i)
    {
        const Vector<int> neighbors = nodeNeighbors[i];

        Real accumArea = area[i];
        for (int j=0; j<neighbors.size(); ++j)
            accumArea += area[ neighbors[j] ];

        Real accumWt = vals[i] * area[i];
        for (int j=0; j<neighbors.size(); ++j)
            accumWt += vals[ neighbors[j] ] * area[ neighbors[j] ];

        newVals[i] = accumWt / accumArea;
    }
}

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    if (ParallelDescriptor::NProcs()>1)
        Abort("Code is not yet parallel safe");

    ParmParse pp;

    std::string infile; pp.get("infile",infile);

    int verbose=0; pp.query("verbose",verbose);

    int Nlev = ml_Nlevels(infile);
    Vector<MultiFab*> streamlines(Nlev);
    Vector<std::string> names;

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cerr << "Reading stream file...\n";

    std::map<std::string,int> idxOfStrInMem;

    int nComp; pp.countval("comps",nComp);
    Vector<int> comps(nComp); pp.getarr("comps",comps);

    //const int nComp(5); // HOW CAN THIS POSSIBLY BE HARDWIRED?

    Vector<int> strComps(nComp);
    for (int i=0; i<nComp; i++)
      strComps[i]=comps[i]+BL_SPACEDIM;

    read_ml_streamline_names(infile,strComps,names);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cerr << "...will read the following components: ";

        for (int i=0; i<names.size(); ++i)
            std::cerr << names[i] << ' ';
        std::cerr << '\n';
    }

    // Find the coordinates in the components
    Vector<int> idX(BL_SPACEDIM,-1);
    for (int i=0; i<names.size(); ++i)
    {
        if (names[i]==std::string("X")) idX[0]=i;
        if (names[i]==std::string("Y")) idX[1]=i;
#if BL_SPACEDIM==3
        if (names[i]==std::string("Z")) idX[2]=i;
#endif
    }

    // Build list of output variables
    int nOutComp=nComp+BL_SPACEDIM;
    std::vector<std::string> outNames(nOutComp);
    std::cerr << "outNames: ";
    for (int i=0; i<nOutComp; ++i) {
        outNames[i] = names[i];
        std::cerr << outNames[i] << " ";
    }
    std::cerr << std::endl;

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cerr << "Reading stream file data: " << infile << "...\n";            

    Vector<int> faceData;
    int nElts;
    Vector<Vector<Vector<int> > > inside_nodes;
    read_ml_streamline_data(infile,strComps,streamlines,names,faceData,nElts,inside_nodes);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cerr << "...finished reading stream file \n";
        std::cerr << "   got the following components: ";

        for (int i=0; i<names.size(); ++i)
            std::cerr << names[i] << ' ';
        std::cerr << '\n';

        std::cerr << "nElts: " << nElts << '\n';
    }

    const int nodesPerElt = faceData.size() / nElts;
    BL_ASSERT(faceData.size()%nodesPerElt==0);

    // Build a structure that for each node points to where in the streamlines to get the data
    std::vector<MLloc> nodeMap = build_nodeMap(inside_nodes);

    // Write a dat file for each variable
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cerr << "Writing output...\n";

    int nComps=outNames.size();
    int nNodes=nodeMap.size();
    // loop over each variable
    for (int iComp=0; iComp<nComps; iComp++) {
      std::cerr << names[iComp] << std::endl;
      std::string outFileName=infile+"/"+names[iComp]+".dat";
      FILE *outfile=fopen(outFileName.c_str(),"w");
      // loop over each node
      for (int iNode=0; iNode<nNodes; iNode++) {
	const MLloc& p=nodeMap[iNode];
	const FArrayBox& fab = (*streamlines[p.amr_lev])[p.box_idx];
	const int nPts = fab.box().length(1);
	IntVect loc = IntVect(D_DECL(p.pt_idx,fab.box().smallEnd()[1],0));
	// loop over each point
	for (int iPt=0; iPt<nPts; iPt++, loc += amrex::BASISV(1)) {
	  Real val = fab(loc,iComp);
	  fprintf(outfile,"%e ",val);
	}
	fprintf(outfile,"\n");
      }
      fclose(outfile);
    }

    amrex::Finalize();
    return 0;
}

int
get_jlo(const Vector<MultiFab*>& data)
{
    int jlo = (*data[0])[0].box().smallEnd(1);
    for (int i=0; i<data.size(); ++i)
        for (int j=0; j<data[i]->size(); ++j)
            jlo = std::min(jlo,(*data[i])[j].box().smallEnd(1));
    return jlo;
}

int
get_nPts(const Vector<MultiFab*>& data)
{
    int nPts = 0;
    for (int i=0; i<data.size(); ++i)
        for (int j=0; j<data[i]->size(); ++j)
            nPts = std::max(nPts,(*data[i])[j].box().length(1));
    return nPts;
}

Real
tetVol(const Real* A,
       const Real* B,
       const Real* C,
       const Real* D)
{
    Real R1[3], R2[3], R3[3], R4[3];
    for (int i=0; i<3; ++i)
    {
        R1[i] = D[i] - A[i];
        R2[i] = B[i] - A[i];
        R3[i] = C[i] - A[i];
    }
    for (int i=0; i<3; ++i)
    {
        R4[0] = R2[1]*R3[2] - R3[1]*R2[2];
        R4[1] = R2[2]*R3[0] - R3[2]*R2[0];
        R4[2] = R2[0]*R3[1] - R3[0]*R2[1];
    }
    Real res=0;
    for (int i=0; i<3; ++i)
        res += R1[i]*R4[i];
    return std::abs(res);
}

Real
max_grad(const MLloc&            p,
         int                     comp,
         const Vector<MultiFab*>& data,
         const Vector<int>&       idX)
{
    const FArrayBox& fab = (*data[p.amr_lev])[p.box_idx];    
    BL_ASSERT(comp<fab.nComp());

    const int nPts = fab.box().length(1);
    Real gradMax = 0.;

    IntVect hi, lo;
    Real hiVal, loVal, hiX[3], loX[3], tot;

    hi = IntVect(D_DECL(p.pt_idx,fab.box().smallEnd()[1],0));
    BL_ASSERT(fab.box().contains(hi));
    hiVal = fab(hi,comp);
    for (int j=0; j<BL_SPACEDIM; ++j)
        hiX[j] = fab(hi,idX[j]);

    // Find an eps
    Real eps, maxs = 0;
    for (int i=1; i<nPts; ++i)
    {
        lo = hi;
        for (int j=0; j<BL_SPACEDIM; ++j)
            loX[j] = hiX[j];

        hi = lo + amrex::BASISV(1);
        BL_ASSERT(fab.box().contains(hi));

        tot = 0;
        for (int j=0; j<BL_SPACEDIM; ++j)
        {
            hiX[j] = fab(hi,idX[j]);
            Real dx = hiX[j] - loX[j];
            tot += dx*dx;
        }
        const Real L = std::sqrt(tot);
        
        maxs = (i==1 ? L : std::max(maxs,L));
    }
    eps = 1.e-4*maxs;

    hi = IntVect(D_DECL(p.pt_idx,fab.box().smallEnd()[1],0));
    int loc = 1;
    for (int i=1; i<nPts; ++i)
    {
        lo = hi;
        loVal = hiVal;
        for (int j=0; j<BL_SPACEDIM; ++j)
            loX[j] = hiX[j];

        hi = lo + amrex::BASISV(1);
        BL_ASSERT(fab.box().contains(hi));
        hiVal = fab(hi,comp);
        tot = 0;
        for (int j=0; j<BL_SPACEDIM; ++j)
        {
            hiX[j] = fab(hi,idX[j]);
            Real dx = hiX[j] - loX[j];
            tot += dx*dx;
        }
        const Real L = std::sqrt(tot);

        if (L>maxs)
        {
            Real grad = std::abs( (hiVal - loVal)/L );
            if (grad >= gradMax) {
                gradMax = grad;
                loc = i;
            }
        }       
    }

    return gradMax;
}

void
peak_val(const MLloc&            p,
         int                     pComp,
         Vector<int>&             sampleComps,
         Vector<Real>&            peakSamples,
         bool*                   peakOK,
         const Vector<MultiFab*>& data)
{
    const FArrayBox& fab = (*data[p.amr_lev])[p.box_idx];    
    BL_ASSERT(pComp<fab.nComp());
    for (int i=0; i<sampleComps.size(); ++i) {
        BL_ASSERT(sampleComps[i]<fab.nComp());
    }
    BL_ASSERT(sampleComps.size()==peakSamples.size());
    const int nPts = fab.box().length(1);

    int peakValLoc = fab.box().smallEnd()[1];
    IntVect iv = IntVect(D_DECL(p.pt_idx,peakValLoc,0));
    BL_ASSERT(fab.box().contains(iv));
    Real peakVal = fab(iv,pComp);
    for (int i=1; i<nPts; ++i)
    {
        iv += amrex::BASISV(1);
        BL_ASSERT(fab.box().contains(iv));
        const Real newVal = fab(iv,pComp);
        if (newVal > peakVal)
        {
            peakVal = newVal;
            peakValLoc = iv[1];
        }
    }

    iv[1] = peakValLoc;
    for (int i=0; i<sampleComps.size(); ++i) {
        peakSamples[i] = fab(iv,sampleComps[i]);
    }

    if ( (peakValLoc == fab.box().smallEnd()[1])
         || (peakValLoc == fab.box().bigEnd()[1]) )
    {
        std::cerr << "peakVal on end of line!" << std::endl;
        *peakOK = false;
    }
    else
    {
        *peakOK = true;
    }
}

Real
wedge_volume_int(const Vector<const MLloc*>& p,
                 int                        ptOnStr,
                 int                        comp,
                 const Vector<MultiFab*>&    data,
                 const Vector<int>&          idX)
{
    BL_ASSERT(comp<data[0]->nComp());

    Real result;

    const int L1 = p[0]->amr_lev;
    const int L2 = p[1]->amr_lev;

    const int B1 = p[0]->box_idx;
    const int B2 = p[1]->box_idx;

    const int P1 = p[0]->pt_idx;
    const int P2 = p[1]->pt_idx;

    if (p.size() == 2)
    {
        Real A[2], B[2], C[2], D[2];

        const IntVect iv1(D_DECL(P1,ptOnStr,0));
        const IntVect iv2(D_DECL(P2,ptOnStr,0));
        const IntVect iv1p(D_DECL(P1,ptOnStr+1,0));
        const IntVect iv2p(D_DECL(P2,ptOnStr+1,0));

        A[0] = (*data[L1])[B1](iv1,idX[0]);
        A[1] = (*data[L1])[B1](iv1,idX[1]);

        B[0] = (*data[L2])[B2](iv2,idX[0]);
        B[1] = (*data[L2])[B2](iv2,idX[1]);

        C[0] = (*data[L1])[B1](iv1p,idX[0]);
        C[1] = (*data[L1])[B1](iv1p,idX[1]);

        D[0] = (*data[L2])[B2](iv2p,idX[0]);
        D[1] = (*data[L2])[B2](iv2p,idX[1]);

        const Real area1 = std::abs(.5 * (((B[0]-A[0])*(D[1]-A[1])) - ((B[1]-A[1])*(D[0]-A[0]))));

        const Real area2 = std::abs(.5 * (((D[0]-A[0])*(C[1]-A[1])) - ((D[1]-A[1])*(C[0]-A[0]))));

        if (comp<0)
        {
            result = area1 + area2;
        }
        else
        {
            Real val1 = 
                (*data[L1])[B1](iv1,comp) + 
                (*data[L2])[B2](iv2,comp) +
                (*data[L2])[B2](iv2p,comp);

            Real val2 = 
                (*data[L1])[B1](iv1,comp) + 
                (*data[L2])[B2](iv2p,comp) +
                (*data[L1])[B1](iv1p,comp);

            val1 /= 3;
            val2 /= 3;

            result = area1*val1 + area2*val2;
            
        }
    }
    else
    {
        BL_ASSERT(p.size() == 3);

        int L3 = p[2]->amr_lev;
        int B3 = p[2]->box_idx;
        int P3 = p[2]->pt_idx;

        Real A[3], B[3], C[3], D[3], E[3], F[3];

        const IntVect iv1(D_DECL(P1,ptOnStr,0));
        const IntVect iv2(D_DECL(P2,ptOnStr,0));
        const IntVect iv3(D_DECL(P3,ptOnStr,0));
        const IntVect iv1p(D_DECL(P1,ptOnStr+1,0));
        const IntVect iv2p(D_DECL(P2,ptOnStr+1,0));
        const IntVect iv3p(D_DECL(P3,ptOnStr+1,0));

        A[0] = (*data[L1])[B1](iv1,idX[0]);
        A[1] = (*data[L1])[B1](iv1,idX[1]);
        A[2] = (*data[L1])[B1](iv1,idX[2]);

        B[0] = (*data[L2])[B2](iv2,idX[0]);
        B[1] = (*data[L2])[B2](iv2,idX[1]);
        B[2] = (*data[L2])[B2](iv2,idX[2]);

        C[0] = (*data[L3])[B3](iv3,idX[0]);
        C[1] = (*data[L3])[B3](iv3,idX[1]);
        C[2] = (*data[L3])[B3](iv3,idX[2]);

        D[0] = (*data[L1])[B1](iv1p,idX[0]);
        D[1] = (*data[L1])[B1](iv1p,idX[1]);
        D[2] = (*data[L1])[B1](iv1p,idX[2]);

        E[0] = (*data[L2])[B2](iv2p,idX[0]);
        E[1] = (*data[L2])[B2](iv2p,idX[1]);
        E[2] = (*data[L2])[B2](iv2p,idX[2]);

        F[0] = (*data[L3])[B3](iv3p,idX[0]);
        F[1] = (*data[L3])[B3](iv3p,idX[1]);
        F[2] = (*data[L3])[B3](iv3p,idX[2]);

        // The following are actually 6 times the tet volume
        const Real vol_EABC = tetVol(A,B,C,E);
        const Real vol_ADEF = tetVol(A,D,E,F);
        const Real vol_ACEF = tetVol(C,E,F,A);

        if (comp<0)
        {
            result = (vol_EABC + vol_ADEF + vol_ACEF)/6.;
        }
        else
        {
            // The following are actually 6 times the tet volume
            const Real vol_DABC = tetVol(A,B,C,D);
            const Real vol_FABC = tetVol(A,B,C,F);
            const Real vol_BDEF = tetVol(B,D,E,F);
            const Real vol_CDEF = tetVol(C,D,E,F);
            const Real vol_ACED = tetVol(C,E,D,A);
            const Real vol_BCDF = tetVol(B,C,D,F);
            const Real vol_BCDE = tetVol(B,C,D,E);
            const Real vol_ABDF = tetVol(B,D,F,A);
            const Real vol_ABEF = tetVol(B,E,F,A);

            const Real vA = (*data[L1])[B1](IntVect(D_DECL(P1,ptOnStr,  0)),comp);
            const Real vB = (*data[L2])[B2](IntVect(D_DECL(P2,ptOnStr,  0)),comp);
            const Real vC = (*data[L3])[B3](IntVect(D_DECL(P3,ptOnStr,  0)),comp);                    
            const Real vD = (*data[L1])[B1](IntVect(D_DECL(P1,ptOnStr+1,0)),comp);
            const Real vE = (*data[L2])[B2](IntVect(D_DECL(P2,ptOnStr+1,0)),comp);
            const Real vF = (*data[L3])[B3](IntVect(D_DECL(P3,ptOnStr+1,0)),comp);

            // These are actually 24 times the integral
            const Real int_1 = ( (vD+vA+vB+vC)*vol_DABC + 
                                 (vB+vD+vE+vF)*vol_BDEF + 
                                 (vB+vC+vD+vF)*vol_BCDF );

            const Real int_2 = ( (vD+vA+vB+vC)*vol_DABC + 
                                 (vC+vD+vE+vF)*vol_CDEF + 
                                 (vB+vC+vD+vE)*vol_BCDE );

            const Real int_3 = ( (vE+vA+vB+vC)*vol_EABC + 
                                 (vA+vD+vE+vF)*vol_ADEF + 
                                 (vA+vC+vE+vF)*vol_ACEF );

            const Real int_4 = ( (vE+vA+vB+vC)*vol_EABC + 
                                 (vC+vD+vE+vF)*vol_CDEF + 
                                 (vA+vC+vE+vD)*vol_ACED );

            const Real int_5 = ( (vF+vA+vB+vC)*vol_FABC + 
                                 (vA+vD+vE+vF)*vol_ADEF + 
                                 (vA+vB+vE+vF)*vol_ABEF );

            const Real int_6 = ( (vF+vA+vB+vC)*vol_FABC + 
                                 (vB+vD+vE+vF)*vol_BDEF + 
                                 (vA+vB+vD+vF)*vol_ABDF );

            // integrals need to be scaled by 1/24, but also must scale average by 6
            result = (int_1 + int_2 + int_3 + int_4 + int_5 + int_6)/144.;
        }
    }

    return result;
}

Real
wedge_surf_area(Vector<const MLloc*> p,
                const Vector<MultiFab*>& data,
                const Vector<int>& idX, int ptOnStr)
{
    const int nComp = p.size();

    BL_ASSERT(nComp>1);
    BL_ASSERT(idX.size()>1);
    BL_ASSERT(idX[0]<nComp);
    BL_ASSERT(idX[1]<nComp);

    const int L1 = p[0]->amr_lev;
    const int L2 = p[1]->amr_lev;
    
    const int B1 = p[0]->box_idx;
    const int B2 = p[1]->box_idx;
    
    const int P1 = p[0]->pt_idx;
    const int P2 = p[1]->pt_idx;
        
    Real res=-1;
  
    if (nComp==2)
    {
        Real A[2], B[2];
        const IntVect iv1(D_DECL(P1,ptOnStr,0));
        const IntVect iv2(D_DECL(P2,ptOnStr,0));
        
        A[0] = (*data[L1])[B1](iv1,idX[0]);
        A[1] = (*data[L1])[B1](iv1,idX[1]);
        
        B[0] = (*data[L2])[B2](iv2,idX[0]);
        B[1] = (*data[L2])[B2](iv2,idX[1]);

        res = std::sqrt( (A[0]-B[0])*(A[0]-B[0]) + (A[1]-B[1])*(A[1]-B[1]) );
    }
    else
    {
        BL_ASSERT(idX.size()>2);
        BL_ASSERT(idX[2]<nComp);
        BL_ASSERT(nComp>2);
        int L3 = p[2]->amr_lev;
        int B3 = p[2]->box_idx;
        int P3 = p[2]->pt_idx;
        
        Real A[3], B[3], C[3];
        const IntVect iv1(D_DECL(P1,ptOnStr,0));
        const IntVect iv2(D_DECL(P2,ptOnStr,0));
        const IntVect iv3(D_DECL(P3,ptOnStr,0));
        
        A[0] = (*data[L1])[B1](iv1,idX[0]);
        A[1] = (*data[L1])[B1](iv1,idX[1]);
        A[2] = (*data[L1])[B1](iv1,idX[2]);
        
        B[0] = (*data[L2])[B2](iv2,idX[0]);
        B[1] = (*data[L2])[B2](iv2,idX[1]);
        B[2] = (*data[L2])[B2](iv2,idX[2]);
        
        C[0] = (*data[L3])[B3](iv3,idX[0]);
        C[1] = (*data[L3])[B3](iv3,idX[1]);
        C[2] = (*data[L3])[B3](iv3,idX[2]);
        
        Real R1[3], R2[3], R3[3];
        for (int i=0; i<3; ++i)
        {
            R1[i] = B[i] - A[i];
            R2[i] = C[i] - A[i];
        }
        R3[0] = R1[1]*R2[2] - R2[1]*R1[2];
        R3[1] = R1[2]*R2[0] - R2[2]*R1[0];
        R3[2] = R1[0]*R2[1] - R2[0]*R1[1];
        res=0;
        for (int i=0; i<3; ++i)
        {
            res += R3[i]*R3[i];
        }
        res = 0.5*std::sqrt(res);
    }
    return res;
}

std::vector<MLloc>
build_nodeMap(const Vector<Vector<Vector<int> > >& inside_nodes)
{
    // Count number of nodes total
    int num_nodes = 0;
    int Nlev = inside_nodes.size();
    for (int i=0; i<Nlev; ++i)
        for (int j=0; j<inside_nodes[i].size(); ++j)
            num_nodes += inside_nodes[i][j].size();

    std::vector<MLloc> nodeMap(num_nodes);
    for (int i=0; i<Nlev; ++i)
    {
        for (int j=0; j<inside_nodes[i].size(); ++j)
        {
            for (int k=0; k<inside_nodes[i][j].size(); ++k)
            {
                BL_ASSERT(inside_nodes[i][j][k]<=num_nodes);
                BL_ASSERT(inside_nodes[i][j][k]>0);
                nodeMap[inside_nodes[i][j][k] - 1] = MLloc(i,j,k); // Remember that inside_nodes is 1-based
            }
        }
    }

    return nodeMap;    
}

int
ml_Nlevels(const std::string& infile)
{
    // Read header for path traces
    const string Header = infile + std::string("/Header");
    std::ifstream ifs;
    ifs.open(Header.c_str());
    std::string typeLabel; ifs >> typeLabel;
    int NlevPath; ifs >> NlevPath;
    ifs.close();
    return NlevPath;
}

static
void
ReadMF (MultiFab&    mf,
        Vector<int>&  comps,
        std::string& infile)
{
    VisMF vismf(infile);

    BL_ASSERT(comps.size()<=vismf.nComp());

    const BoxArray& ba = vismf.boxArray();
    const DistributionMapping dm(ba);
    mf.define(ba,dm,comps.size(),vismf.nGrow());

    for (int j=0; j<comps.size(); ++j)
    {
        BL_ASSERT(comps[j]>=0 && comps[j]<vismf.nComp());
    }

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        for (int j=0; j<comps.size(); ++j)
        {
            const FArrayBox& fab = vismf.GetFab(mfi.index(), comps[j]);
            mf[mfi].copy(fab,0,j,1);
            vismf.clear(mfi.index(), comps[j]);
        }
    }
    BL_ASSERT(mf.ok());
}

void
read_ml_streamline_names(const std::string& infile,
                         Vector<int>&        comps,
                         vector<string>&    names)
{
    // Read header for path traces
    const string Header = infile + string("/Header");
    std::ifstream ifh;
    ifh.open(Header.c_str(), std::ifstream::in);
    if (!ifh.is_open()) {
      std::cout << " *** Error opening Header file *** " << Header << std::endl;
      exit(1);
    }
    std::string typeLabel; ifh >> typeLabel;
    int NlevPath; ifh >> NlevPath;
    int nCompPath; ifh >> nCompPath;
    vector<string> fileNames(nCompPath);
    for (int i=0; i<nCompPath; ++i)
        ifh >> fileNames[i];
    ifh.close();

    std::cout << "NlevPath:  " << NlevPath  << '\n';
    std::cout << "nCompPath: " << nCompPath << '\n';

    // Find the coordinates in the components
    Vector<int> readComps(comps.size()+BL_SPACEDIM,-1);
    names.resize(readComps.size());
    names[0] = "X"; names[1] = "Y";
#if BL_SPACEDIM==3
 names[2] = "Z";
#endif
    for (int i=0; i<nCompPath; ++i)
        for (int j=0; j<BL_SPACEDIM; ++j)
            if (fileNames[i]==names[j]) readComps[j] = i;

    std::cout << " nCompPath " << nCompPath << std::endl;
    for (int i=0; i<nCompPath; ++i)
      std::cout << " fileNames[i] " << fileNames[i] << std::endl;

    std::cout << " comps.size() = " << comps.size() << std::endl;
    for (int i=0; i<comps.size(); ++i)
      std::cout << comps[i] << std::endl;
    
    for (int i=0; i<comps.size(); ++i)
    {
        std::cout << " i = " << i << std::endl;
        readComps[BL_SPACEDIM+i] = comps[i];
        std::cout << " readComps[BL_SPACEDIM+i] = " << readComps[BL_SPACEDIM+i] << std::endl;
        names[BL_SPACEDIM+i] = fileNames[comps[i]];
        std::cout << " names[BL_SPACEDIM+i] = " << names[BL_SPACEDIM+i] << std::endl;
    }

}

void
read_ml_streamline_data(const std::string& infile,
                        Vector<int>&        comps,
                        Vector<MultiFab*>&  data,
                        vector<string>&    names,
                        Vector<int>&        faceData,
                        int&               nElts,
                        Vector<Vector<Vector<int> > >& inside_nodes)
{
    // Read header for path traces
    const string Header = infile + string("/Header");
    std::ifstream ifh;
    ifh.open(Header.c_str());
    std::string typeLabel; ifh >> typeLabel;
    int NlevPath; ifh >> NlevPath;
    int nCompPath; ifh >> nCompPath;
    vector<string> fileNames(nCompPath);
    for (int i=0; i<nCompPath; ++i)
        ifh >> fileNames[i];
    ifh.close();

    std::cout << "NlevPath:  " << NlevPath  << '\n';
    std::cout << "nCompPath: " << nCompPath << '\n';

    // Find the coordinates in the components
    Vector<int> readComps(comps.size()+BL_SPACEDIM,-1);
    names.resize(readComps.size());
    names[0] = "X"; names[1] = "Y";
#if BL_SPACEDIM==3
 names[2] = "Z";
#endif
    for (int i=0; i<nCompPath; ++i)
        for (int j=0; j<BL_SPACEDIM; ++j)
            if (fileNames[i]==names[j]) readComps[j] = i;
    for (int i=0; i<comps.size(); ++i)
    {
        readComps[BL_SPACEDIM+i] = comps[i];
        names[BL_SPACEDIM+i] = fileNames[comps[i]];
    }

    // Read elements
    const string iElements = infile + string("/Elements");
    std::ifstream ife;
    ife.open(iElements.c_str());
    ife >> nElts;
    int nodesPerElt; ife >> nodesPerElt;
    faceData.resize(nElts*nodesPerElt);
    for (int i=0; i<nElts*nodesPerElt; ++i)
        ife >> faceData[i];

    for (int lev=0; lev<NlevPath; ++lev)
    {
        char buf[64];
        sprintf(buf, "/Level_%d", lev);
        std::string ThisStreamFile = infile + std::string(buf);
        ThisStreamFile += "/Str";
        data[lev] = new MultiFab();
        std::cout << "Calling ReadMF() at lev: " << lev << " ...\n";
        ReadMF(*data[lev],readComps,ThisStreamFile);
    }
    //
    // Force other processors to wait until parallel data has been read
    //
    ParallelDescriptor::Barrier();

    // Read element distribution (after we know Nlev from path traces)
    inside_nodes.resize(NlevPath);
    for (int lev=0; lev<NlevPath; ++lev)
    {
        inside_nodes[lev].resize(data[lev]->size());

        // Get number of nonempty boxes
        int num_non_zero;
        ife >> num_non_zero;

        for (int j=0; j<num_non_zero; ++j)
        {
            // Get index of non-empty box, and number of entries
            int box_id, num_ids;
            ife >> box_id >> num_ids;

            inside_nodes[lev][box_id].resize(num_ids);

            for (int k=0; k<num_ids; ++k)
                ife >> inside_nodes[lev][box_id][k];
        }
    }
    ife.close();
}

#ifdef BL_WRITE_BINARY
void
write_binary_tec_file(const std::string&          outfile,
                      const std::vector<string>&  names,
                      const Vector<int>&           faceData,
                      const Vector<MultiFab*>&     nodes,
                      const std::vector<MLloc>&   nodeMap,
                      const Vector<Vector<Real> >&  eltData,
                      const std::string&          title,
                      const Vector<int>&           idX)
{                
    std::string vars = "X Y Z";
    for (int j=0; j<names.size(); ++j)
        vars += " " + names[j];

    INTEGER4 I;
    INTEGER4 Debug = 0;
    INTEGER4 VIsDouble = 1;
    I = TECINI100("ELT STATS (simulating element-centered data",
                  (char*)vars.c_str(),
                  (char*)outfile.c_str(),
                  ".",
                  &Debug,
                  &VIsDouble);
    
    INTEGER4 ZoneType  = 2;  /* FETRIANGLE */
    INTEGER4 KMax      = 0;  /* Unused */
    INTEGER4 ICellMax  = 0;
    INTEGER4 JCellMax  = 0;
    INTEGER4 KCellMax  = 0;
    INTEGER4 IsBlock   = 1;
    INTEGER4 NumFaceConnections = 0;
    INTEGER4 FaceNeighborMode   = 0;
    INTEGER4 ShareConnectivityFromZone = 0;

    int nComp = eltData[0].size() + BL_SPACEDIM;
    int nPtsOnStr = get_nPts(nodes);
    int nElts = eltData.size();
    int nodesPerElt = faceData.size() / nElts;
    int nPts = nElts*nodesPerElt;
    
    Box nbox(IntVect::TheZeroVector(),IntVect(nPts-1,0,0));
    FArrayBox nodeFab(nbox,nComp);
    Vector<int> connData(nPts);
    // Build nodeFab so we can write it, use ptOnStr=0 to get point at surface
    // however, do it so that the nodes are multiply defined...my way to show elt-based
    // data (as opposed to node-based data)

    cout << "Building new node data" << endl;
    int ptOnStr = 0;
    for (int i=0; i<nElts; ++i)
    {
        int offset = i*nodesPerElt;
        for (int j=0; j<nodesPerElt; ++j)
        {
            const MLloc& n = nodeMap[ faceData[offset+j] - 1 ];
            const IntVect iv = IntVect(n.pt_idx,ptOnStr,0);
            nodeFab(IntVect(D_DECL(offset+j,0,0)),0) = (*nodes[n.amr_lev])[n.box_idx](iv,idX[0]);
            nodeFab(IntVect(D_DECL(offset+j,0,0)),1) = (*nodes[n.amr_lev])[n.box_idx](iv,idX[1]);
            nodeFab(IntVect(D_DECL(offset+j,0,0)),2) = (*nodes[n.amr_lev])[n.box_idx](iv,idX[2]);
            for (int m=0; m<eltData[0].size(); ++m)
                nodeFab(IntVect(D_DECL(offset+j,0,0)),m+BL_SPACEDIM) = eltData[i][m];
            connData[offset+j] = offset+ j + 1;
        }
    }

    I = TECZNE100((char*)title.c_str(),
                  &ZoneType,
                  &nPts,
                  &nElts,
                  &KMax,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &IsBlock,
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  NULL,      /* ValueLocation */
                  NULL,      /* ShareVarFromZone */
                  &ShareConnectivityFromZone);
    
    INTEGER4 III = nPts * nComp;
    I = TECDAT100(&III,nodeFab.dataPtr(),&VIsDouble);
    I = TECNOD100(connData.dataPtr());
    I = TECEND100(); 
}
#else
void
write_ascii_tec_file(const std::string&           outfile,
                     const std::vector<string>&   names,
                     const Vector<int>&           faceData,
                     const Vector<MultiFab*>&     nodes,
                     const std::vector<MLloc>&    nodeMap,
                     const Vector<Vector<Real> >& eltData,
                     const std::string&           title,
                     const Vector<int>&           idX)
{                
    std::string vars = "VARIABLES = X Y Z";
    for (int j=0; j<names.size(); ++j)
        vars += " " + names[j];

    // Write ASCII output file
    std::ofstream os(outfile.c_str(),std::ios::out);
    int nComp = eltData[0].size() + BL_SPACEDIM;
    int nElts = eltData.size();
    int nodesPerElt = faceData.size() / nElts;
    int nPts = nElts*nodesPerElt;
    
    Box nbox(IntVect::TheZeroVector(),IntVect(D_DECL(nPts-1,0,0)));
    FArrayBox nodeFab(nbox,nComp);
    Vector<int> connData(nPts);
    // Build nodeFab so we can write it, use ptOnStr=0 to get point at surface
    // however, do it so that the nodes are multiply defined...my way to show elt-based
    // data (as opposed to node-based data)

    cout << "Building new node data" << endl;
    int ptOnStr = 0;
    for (int i=0; i<nElts; ++i)
    {
        int offset = i*nodesPerElt;
        for (int j=0; j<nodesPerElt; ++j)
        {
            const MLloc& n = nodeMap[ faceData[offset+j] - 1 ];
            const IntVect iv = IntVect(D_DECL(n.pt_idx,ptOnStr,0));
            for (int k=0; k<BL_SPACEDIM; ++k)
                nodeFab(IntVect(D_DECL(offset+j,0,0)),k) = (*nodes[n.amr_lev])[n.box_idx](iv,idX[k]);
            for (int m=0; m<eltData[0].size(); ++m)
                nodeFab(IntVect(D_DECL(offset+j,0,0)),m+BL_SPACEDIM) = eltData[i][m];
            connData[offset+j] = offset+ j + 1;
        }
    }

    os << vars << std::endl;
    os << "ZONE T=\"" << title << "\" N=" << nPts << " E=" << nElts
       << " F=FEBLOCK ET=" << (nodesPerElt==2 ? "LINESEG" : "TRIANGLE") << std::endl;

    int valsPerLine = 5;
    for (int k=0; k<nComp; ++k)
    {
        const Real* nodeData = nodeFab.dataPtr(k);
        for (int i=0; i<nPts; ++i)
            os << nodeData[i] << (i%valsPerLine==valsPerLine-1 ? "\n" : " ");
        os << std::endl;
    }
    for (int i=0; i<nElts; ++i)
    {
        int offset = i*nodesPerElt;
        for (int k=0; k<nodesPerElt; ++k)
            os << connData[offset+k] << " ";
        os << std::endl;
    }
    os.close();
}
#endif

void
write_binary_mef_file(const std::string&           outfile,
                      const std::vector<string>&   names,
                      const Vector<int>&           faceData,
                      const Vector<MultiFab*>&     nodes,
                      const std::vector<MLloc>&    nodeMap,
                      const Vector<Vector<Real> >& eltData,
                      const std::string&           title,
                      const Vector<int>&           idX)
{                
    std::string vars = "X Y";

#if BL_SPACEDIM==3
    vars += " Z";
#endif

   
    for (int j=0; j<names.size(); ++j)
        vars += " " + names[j];

    int nComp = eltData[0].size() + BL_SPACEDIM;
    int nElts = eltData.size();
    int nodesPerElt = faceData.size() / nElts;


#undef FAKE_NODE_MEF
#define FAKE_NODE_MEF
#ifdef FAKE_NODE_MEF
    int nPts = nElts*nodesPerElt;
    
    Box nbox(IntVect::TheZeroVector(),IntVect(D_DECL(nPts-1,0,0)));
    FArrayBox fab(nbox,nComp);
    Vector<int> connData(nPts);
    // Build fab so we can write it, use ptOnStr=0 to get point at surface
    // however, do it so that the nodes are multiply defined...my way to show elt-based
    // data (as opposed to node-based data)

    cout << "Building new node data" << endl;
    int ptOnStr = 0;
    Real* floatDat = fab.dataPtr();
    for (int i=0; i<nElts; ++i)
    {
        int offset = i*nodesPerElt;
        for (int j=0; j<nodesPerElt; ++j)
        {
            const MLloc& n = nodeMap[ faceData[offset+j] - 1 ];
            const IntVect iv = IntVect(D_DECL(n.pt_idx,ptOnStr,0));
            
            for (int k=0; k<BL_SPACEDIM; ++k)
                *floatDat++ = (*nodes[n.amr_lev])[n.box_idx](iv,idX[k]);

            for (int m=0; m<eltData[0].size(); ++m)
                *floatDat++ = eltData[i][m];
            connData[offset+j] = offset+ j + 1;
        }
    }
#else
    int nPts = nodeMap.size();    
    Box nbox(IntVect::TheZeroVector(),IntVect(nPts-1,0,0));
    FArrayBox fab(nbox,BL_SPACEDIM); // node positions, fab correct size but fill/write transposed below
    Vector<int> connData(nPts);
    int ptOnStr = 0;
    Real* floatDat = fab.dataPtr();
    for (int i=0; i<nPts; ++i)
    {
        const MLloc& n = nodeMap[i];
        const IntVect iv = IntVect(n.pt_idx,ptOnStr,0);            
        for (int k=0; k<BL_SPACEDIM; ++k)
            *floatDat++ = (*nodes[n.amr_lev])[n.box_idx](iv,idX[k]);
    }
    vars += " __data_is_elt_centered";
#endif

    std::ofstream ofs;
    ofs.open(outfile.c_str());

    ofs << title << endl;
    ofs << vars << std::endl;
    ofs << nElts << " " << nodesPerElt << endl;
    fab.writeOn(ofs);
    ofs.write((char*)connData.dataPtr(),sizeof(int)*connData.size());
        
#ifndef FAKE_NODE_MEF
    int nCompElt = eltData[0].size();
    Box ebox(IntVect::TheZeroVector(),IntVect(D_DECL(nElts-1,0,0)));
    fab.resize(ebox,nCompElt);
    floatDat = fab.dataPtr();
    for (int i=0; i<nElts; ++i)
        for (int j=0; j<nCompElt; ++j)
            *floatDat++ = eltData[i][j];
    fab.writeOn(ofs);
#endif

    ofs.close();
}

static
std::vector<std::string> parseVarNames(std::istream& is)
{
    std::string line;
    std::getline(is,line);
    return Tokenize(line,std::string(", "));
}

static std::string parseTitle(std::istream& is)
{
    std::string line;
    std::getline(is,line);
    return line;
}

void
read_iso(const std::string& infile,
         FArrayBox&         nodes,
         Vector<int>&       faceData,
         int&               nElts,
         vector<string>&    names,
         string&            label)
{
    std::ifstream ifs;
    ifs.open(infile.c_str(),std::ios::in|std::ios::binary);
    label = parseTitle(ifs);
    names = parseVarNames(ifs);
    const int nCompSurf = names.size();

    int nodesPerElt;
    ifs >> nElts;
    ifs >> nodesPerElt;

    FArrayBox tnodes;
    tnodes.readFrom(ifs);
    const int nNodes = tnodes.box().numPts();

    // "rotate" the data so that the components are 'in the right spot for fab data'
    nodes.resize(tnodes.box(),nCompSurf);
    Real** np = new Real*[nCompSurf];
    for (int j=0; j<nCompSurf; ++j)
        np[j] = nodes.dataPtr(j);

    Real* ndat = tnodes.dataPtr();
    for (int i=0; i<nNodes; ++i)
    {
        for (int j=0; j<nCompSurf; ++j)
        {
            np[j][i] = ndat[j];
        }
        ndat += nCompSurf;
    }
    delete [] np;
    tnodes.clear();

    faceData.resize(nElts*nodesPerElt,0);
    ifs.read((char*)faceData.dataPtr(),sizeof(int)*faceData.size());
}
