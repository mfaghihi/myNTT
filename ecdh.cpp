#include <NTL/BasicThreadPool.h>
#include <NTL/tools.h>
#include <NTL/vector.h>
#include <NTL/vec_lzz_p.h>
#include <NTL/ZZ.h>
#include <NTL/lzz_p.h>
#include <NTL/lzz_pX.h>
#include <NTL/FFT_impl.h>
#include <NTL/lzz_p.h>

#include <algorithm>
#include <cstdint>
#include <forward_list>
#include <fstream>
#include <iostream>
#include <ostream>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <NTL/ZZ_pXFactoring.h>

using namespace std;
using namespace NTL;

using NTL::ZZ;
using NTL::zz_p;
typedef NTL::zz_p lzz_p;
using NTL::zz_pX;
typedef NTL::zz_pX lzz_pX;
using NTL::vec_zz_p;
typedef NTL::vec_zz_p vec_lzz_p;
typedef NTL::ZZ word;
using NTL::Vec;

//Generates a random number d in the field s.t 0 < d < n-1
ZZ_p KeyGen(ZZ n){
    ZZ_p randZZ_p;
    //ZZ randZZ = RandomBnd(n);
    conv(randZZ_p, RandomBnd(n-1)+1);
    //cout <<n-randZZ<<endl;
    return randZZ_p;
}
// Adds two pints P and Q on the curve ///TODO: It does not check if the points are on the curve
Vec<ZZ_p> pointsAdd(Vec<ZZ_p>P, Vec<ZZ_p>Q){
    Vec<ZZ_p> result;
    result.SetLength(2);
    if (P[0]==0 && P[1]==0){
        return Q;
    }
    if (Q[0]==0 && Q[1]==0){
        return P;
    }
    //P = -Q;
    if (P[0]==Q[0] && P[1]!=Q[1]){
        conv(result[0],"0");
        conv(result[1],"0");
        return result;
    }
    //P = Q
    if (P[0]==Q[0]){
        //cout << "here pointsAdd 1: P[1]= " << P[1] <<endl;
        ZZ_p s= (3*P[0]*P[0] -3)/(2*P[1]);
        result[0] = s*s-2*P[0];
        result[1] = s*(P[0]-result[0])-P[1]; 
        return result;
    }
    //P != Q
    else {
        ZZ_p s = (P[1]-Q[1])/(P[0]-Q[0]);
        result[0] =s*s - (P[0]+Q[0]);
        result[1] = s*(P[0]-result[0])-P[1];

    }
    return result;
}

// Multiplies the pint G on the curve with scalar du  ///TODO: It does not check if the point is on the curve
Vec<ZZ_p>  scalarMult(Vec<ZZ_p>G, ZZ_p du, ZZ_p n){
    Vec<ZZ_p> result;
    result.SetLength(2);
    Vec<ZZ_p> addend = G;

    ZZ ZZdu, ZZn;
    conv(ZZdu,du);
    conv(ZZn,n);

    conv(result[0],"0");
    conv(result[1],"0");
    //cout << "here scalarMult 1"<<endl;

    if (ZZdu % ZZn == 0){
        conv(result[0],"0");
        conv(result[1],"0");
        return result;
    }
    while (ZZdu != 0){

        if (IsOdd(ZZdu)){
        // add
            result = pointsAdd(result,addend);
        }
        // Double
        addend = pointsAdd(addend,addend);

        ZZdu>>=1;

    }

    return result;
}
int main(){
    ZZ thePrime;
    conv(thePrime,"115792089210356248762697446949407573530086143415290314195533631308867097853951");
    ZZ_p::init(thePrime);
    ZZ_p a(-3);
    ZZ_p b;
    conv(b,"41058363725152142129326129780047268409114441015993725554835256314039467401291");
    Vec<ZZ_p> G;
    G.SetLength(2);
    conv(G[0],"48439561293906451759052585252797914202762949526041747995844080717082404635286");
    conv(G[1],"36134250956749795798585127919587881956611106672985015071877198253568414405109");
    ZZ_p n;
    conv(n,"1157920892103562487626974469494075735299969552241357603424222590610685120443");
    ZZ_p h(1);
    
    ZZ_p sampledu1,sampledu2;
    //conv(sampledu,"968010264577893948860466785904642298477494770782054870952873973004363905723");
    conv(sampledu1,"9041875927920587073444478640532130384204131118514724791156633039651898202668");
    conv(sampledu2,"16487785647343886075367572781924531209336180327842376707869205506526978857988");

    ZZ aa;
    conv(aa,sampledu1);
    //cout<< G <<endl;
    cout << "mult result with sampledu1 = "<<scalarMult(G,sampledu1,n)<<endl;
    cout << "mult result with sampledu2 = "<<scalarMult(G,sampledu2,n)<<endl;
    cout << "shared key = "<<scalarMult(scalarMult(G,sampledu1,n),sampledu2,n)<<endl;
    cout << "shared key = "<<scalarMult(scalarMult(G,sampledu2,n),sampledu1,n)<<endl;

    cout<<"Finish"<<endl;
}
