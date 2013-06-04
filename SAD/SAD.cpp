#include <cstdio>
#include <memory>
#include "SAD.h"

int ADS::nvar = 0;
int ADS::nnz_pd = 0;
int ADS::maxNNZ_pd = 0;
int* ADS::cooRow = 0;
int* ADS::cooCol = 0;
float* ADS::pd = 0;

void ADS::Initialize()
{
	Clear();
	nvar = 0;
	nnz_pd = 0;
	maxNNZ_pd = MAX_NVAR;
	cooRow = new int[MAX_NVAR];
	cooCol = new int[MAX_NVAR];
	pd = new float[MAX_NVAR];
}

void ADS::Clear()
{
	if(cooRow)
	{
		delete [] cooRow;
		delete [] cooCol;
		delete [] pd;
	}
	cooRow = 0;
	cooCol = 0;
	pd = 0;
}

void ADS::GetJacobianForward(float *J, int m, int n)
{
	float *adj = new float[ADS::nvar];
	memset(adj, 0, sizeof(*adj)*ADS::nvar);

	int *rid = ADS::cooRow;
	int *cid = ADS::cooCol;
	float *pd = ADS::pd;

	for(int xid = 0; xid < n; ++xid)
	{
		adj[xid] = 1;
		for(int j = 0; j < ADS::nnz_pd; ++j)
		{
			adj[ rid[j] ] += pd[j]*adj[cid[j]];
		}
		memcpy(J+m*xid, adj + ADS::nvar - m, sizeof(*J)*m);
		memset(adj, 0, sizeof(*adj)*ADS::nvar);
	}
	delete [] adj;
}
void ADS::GetJacobianReverse(float *Jt, int m, int n)
{
	float *adj = new float[ADS::nvar];
	memset(adj, 0, sizeof(*adj)*ADS::nvar);

	int *rid = ADS::cooRow;
	int *cid = ADS::cooCol;
	float *pd = ADS::pd;

	for(int yid = 0; yid < m; ++yid)
	{
		int xid = ADS::nvar - m + yid;
		adj[xid] = 1;
		for(int j = ADS::nnz_pd - 1; j >= 0; --j)
		{
			adj[ cid[j] ] += pd[j]*adj[ rid[j] ];
		}
		memcpy(Jt+m*yid, adj, sizeof(*Jt)*n);
		memset(adj, 0, sizeof(*adj)*ADS::nvar);
	}
	delete [] adj;
}

void ADS::ShowNodes()
{
	for(int i = 0; i < nnz_pd; ++i)
	{
		printf("(%2d, %2d)\t%6.2f\n", cooRow[i], cooCol[i], pd[i]);
	}
	putchar('\n');
}

void ShowJacobian(float *J, int m, int n, bool Transposed)
{
	for(int i = 0; i < m; ++i)
	{
		for(int j = 0; j < n; ++j)
		{
			if(Transposed)
			{
				printf("%6.4f\t", J[i*n + j] );
			}
			else{
				printf("%6.4f\t", J[i+j*m] );
			}
		}
		putchar('\n');
	}
}


void TestCase1()
{
	ADS::Initialize();

	ADV U[5], V[5];
	for(int i = 1; i < 4; ++i)
		U[i] = 1;
	U[0] = 1;
	U[4] = -1;

	for(int i = 1; i < 4; ++i)
	{
		V[i] = Operation4(-2, U[i-1], U[i], U[i+1]);
	}

	ADS::ShowNodes();
	putchar('\n');
	printf("%6.4f\t%6.4f\t%6.2f\n",V[1].v, V[2].v, V[3].v);
	putchar('\n');

	float *J = new float[3*3];
	ADS::GetJacobianForward(J, 3, 3);

	ADS::Clear();
}

void TestCase2()
{
	ADS::Initialize();

	int M = 4;
	int N = M*M;
	float L = 1;
	float T = 0.5;

	ADV *s = new ADV[M-1];
	for(int i = 0; i < M-1; ++i)
		s[i] = 1;

	float lambda = T*M*M/(N*L*L);
	ADV *K = new ADV[(M-1)*3], *pK;
	K[0] = Operation1(lambda, s[0], s[1]);
	K[1] = Operation2(lambda, s[0]);
	K[2] = Operation3(lambda, s[0], s[1]);
	for(int i = 1; i < (M-2); ++i)
	{
		pK = K + i*3;
		pK[0] = Operation3(lambda, s[i-1], s[i]);
		pK[1] = Operation4(lambda, s[i-1], s[i], s[i+1]);
		pK[2] = Operation3(lambda, s[i], s[i+1]);
	}
	pK = K + (M-2)*3;
	pK[0] = Operation3(lambda, s[M-3], s[M-2]);
	pK[1] = Operation2(lambda, s[M-2]);
	pK[2] = Operation1(lambda, s[M-2], s[M-3]);

	ADV *ta, *tb, *tc;
	ta = new ADV[M+1];
	tb = new ADV[M+1];
	ta[0] = 1;
	tb[0] = 1;
	for(int i = 1; i < M; ++i)
		ta[i] = 0;
	ta[M] = 0.5;
	tb[M] = 0.5;

	for(int n = 0; n < N; ++n)
	{
		for(int j = 1; j < M; ++j)
		{
			pK = K + (j-1)*3;
			tb[j] = InnerProd3(pK, ta+j-1);
		}
		tc = tb;
		tb = ta;
		ta = tc;
	}

	float u[3] = {0};
	ADV sqe = SquaredError(3, u, ta+1);

	for(int j = 1; j < M; ++j)
	{
		pK = K + (j-1)*3;
		printf("%6.4f\t%6.4f\t%6.4f\n", pK[0].v, pK[1].v, pK[2].v);
	}

	putchar('\n');
	for(int i = 0; i < (M+1); ++i)
	{
		printf("(%d)\t%6.4f\n", ta[i].id, ta[i].v);
	}
	putchar('\n');
	
	printf("sqe = %6.4f\n\n", sqe.v);

	putchar('\n');

	float *J = new float[3*3];
	ADS::GetJacobianForward(J, 1, 3);
	ShowJacobian(J,1,3);
	putchar('\n');
	ADS::GetJacobianReverse(J, 1, 3);
	ShowJacobian(J,1,3,true);
	ADS::Clear();
	delete [] J;
}

