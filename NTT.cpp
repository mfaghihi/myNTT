#include <NTL/BasicThreadPool.h>
#include <NTL/FFT.h>
#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>
#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/tools.h>
#include <NTL/vector.h>
#include <NTL/vec_lzz_p.h>
#include <NTL/ZZ.h>
#include <NTL/lzz_p.h>
#include <NTL/lzz_pX.h>


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

ZZ_p findnthRoot(ZZ_p n, ZZ p){
	ZZ_p result;
	random(result);
	while(result==ZZ_p(0))
	{
		random(result);
	}
	ZZ n1;
    conv(n1,n);
	// If n|(p-1), then there exists n-th root of unity x which must satisfy x^((p-1)/n)!=1.
	// Later we check if x^t != 1 for t < n.
	while (power(result,(p-1)/n1) == ZZ_p(1))
	{
		random(result);
		cout <<result<<endl;
	}

	return power(result,(p-1)/n1);
	
}

bool isPrimRootOfUnity(ZZ_p n,ZZ_p root)
{

	long n1;
	conv(n1,n);
	// Check if x^t != 1 for t < n
	for (long i =1; i<n1;i++)
	{

		if (power(root,i)==1)
		{

			return false;
		}
	}
	cout<< root <<" is "<< n <<"th root of unity"<<endl;
	return true;
}

//void Printres(vector)

int main()
{
   //ZZ p(4611686018326724609);
   ZZ p(4611686018309947393);
   //ZZ p(11);
   //cin >> p;
   ZZ_p::init(p);

   // size of input vector to NTT
   long n1=512; 
   //long n1=5;
	ZZ_p n(n1);
   //ZZ_p Nm1(4611686018326724608);
   ZZ_p Nm1(4611686018309947392); // p-1 (or N-1)
   //ZZ_p Nm1(10);
   //cin >> Nm1;
	ZZ_p k=Nm1/n;
	//conv(k,576460752290840576);
   //Vec< Pair< ZZ_pX, long > > factors;

   //CanZass(factors, f);  // calls "Cantor/Zassenhaus" algorithm

   cout <<"f/"<< n <<"= "<< Nm1/n << "\n";
   cout<<"k= "<<k<<endl;
   ZZ_p nth=findnthRoot(n,p);

   //make sure nth is a primitive root of unity
   while (!isPrimRootOfUnity(n,nth))
   {
	   nth=findnthRoot(n,p);
   }
	
	// Input vector
	Vec<ZZ_p> X;
	X.SetLength(n1);
	/*X[0]= ZZ_p(6);
	X[1]= ZZ_p(0);
	X[2]= ZZ_p(10);
	X[3]= ZZ_p(7);
	X[4]= ZZ_p(2);*/
	
	//nth=ZZ_p(3);

	//Output vector of NTT
	Vec<ZZ_p> Y;
	Y.SetLength(n1);

	//Compute NTT
	ZZ_p tmp;
	for (long i = 0; i<n1; i++)
	{
		X[i] = (random_ZZ_p());
	}
	cout <<"X ="<< X<<endl;
		for (long k = 0; k<n1; k++)
	{
		tmp = 0;
			for (long i = 0; i<n1; i++)
		{
			tmp += X[i]*power(nth,i*k);
		}
		Y[k]=tmp;
	}
	cout <<endl<<"+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+="<<endl;
	cout <<"+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+="<<endl;
	cout <<"+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+="<<endl<<endl;

	cout <<"Y ="<< Y<<endl;

	//Output vector of invNTT
	Vec<ZZ_p> Z;
	Z.SetLength(n1);
	//Compute invNTT
	for (long k = 0; k<n1; k++)
	{
		tmp = 0;
		for (long i = 0; i<n1; i++)
		{
			tmp += Y[i]*power(nth,-i*k);
		}
		Z[k]=tmp/n1;
	}
	//cout <<"Z ="<< Z<<endl;
	
	cout <<"Z-X ="<< Z-X<<endl;
	return 0;
	
}