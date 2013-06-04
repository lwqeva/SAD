#include "HeatConduction.h"
#include "SAD.h"

#include <cmath>
#include <cstdio>
#define M_PI       3.14159265358979323846


void TestCase1()
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

//	float u[3] = {0};
//	ADV sqe = SquaredError(3, u, ta+1);

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
	
//	printf("sqe = %6.4f\n\n", sqe.v);

	putchar('\n');

	float *J = new float[3*3];
	ADS::GetJacobianForward(J, 3, 3);
	ShowJacobian(J,3,3);
	putchar('\n');
	ADS::GetJacobianReverse(J, 3, 3);
	ShowJacobian(J,3,3,true);
	ADS::Clear();
	delete [] J;
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

void TestCase3()
{
	Rod::Initialize(4);

	float *s = new float[Rod::M-1];
	for(int i = Rod::M - 2; i >= 0; --i)
		s[i] = 1;
	Rod rod = Rod::CreateRod(s);
	
	ADS::Initialize();

	ADV *S = new ADV[Rod::M-1];
	for(int i = 0; i < (Rod::M-1); ++i)
	{
		S[i] = s[i];
	}

	ADV *tN = new ADV[Rod::M+1];
	Solve_HeatEquation(S, tN);
	
	for(int i = Rod::M - 2; i >= 0; --i)
		s[i] = 0;
	ADV sqe = SquaredError(Rod::M-1, s, tN+1);

	float *J = new float[(Rod::M-1)*(Rod::M-1)];
	ADS::GetJacobianForward(J, 1, Rod::M-1);
	ShowJacobian(J, 1, Rod::M-1);	

	for(int i = 1; i < Rod::M; ++i)
	{
		rod.tN[i] = tN[i].v;
	}
	rod.Show();

	delete [] s;
	delete [] J;
	delete [] S;
	delete [] tN;
	rod.Destroy();
	ADS::Clear();
	Rod::Clear();
}


int main()
{
	TestCase4();
	return 0;
}
