#include "HeatConduction.h"
#include "SAD.h"
#include <memory>

float Rod::L = 0;
float Rod::T = 0;
int Rod::M = 0;
int Rod::N = 0;
float Rod::lambda = 0;
float* Rod::t0 = 0;

void Rod::Initialize(int M)
{
	L = 1;
	T = 0.5;
	Rod::M = M;
	N = M*M;
	lambda = T*M*M/(N*L*L);
	t0 = new float[M+1];
	memset(t0, 0, sizeof(*t0)*(M+1));
	t0[0] = 1;
	t0[M] = 0.5;
}

void Rod::Clear()
{
	if(t0)
		delete [] t0;
}

Rod Rod::CreateRod(float *s, float *tN)
{
	Rod rod;
	rod.s = new float[Rod::M-1];
	if(s)    // conductivity of interior points
	{
		memcpy(rod.s, s, sizeof(*s)*(Rod::M-1));
	}else
	{
		memset(rod.s, 0, sizeof(*s)*(Rod::M-1));
	}

	rod.tN = new float[Rod::M+1];
	if(tN)
	{
		memcpy(rod.tN, tN, sizeof(*tN)*(Rod::M+1));
	}else
	{
		memset(rod.tN, 0, sizeof(*tN)*(Rod::M+1));
	}
	rod.tN[0] = Rod::t0[0];
	rod.tN[Rod::M] = Rod::t0[M];
	return rod;
}

void Rod::Destroy()
{
	delete [] s;
	delete [] tN;
	s = 0;
	tN = 0;
}

void Rod::Show(bool isTerminal)
{
	float *t = t0;
	if( isTerminal )
		t = tN;
	putchar('\n');
	for(int i = 0; i <= M; ++i)
		printf("%6.4f\n", t[i]);
	putchar('\n');
}

template <class Type>
void Solve_HeatEquation(Type *s, Type *tN)
{
	// Prepare constants and allocate memory for variables
	int M = Rod::M;
	int N = Rod::N;
	float lambda = Rod::lambda;

	Type *K = new Type[(M-1)*3], *pK;
	K[0] = Operation1(lambda, s[0], s[1]);
	K[1] = Operation2(lambda, s[0]);
	K[2] = Operation3(lambda, s[0], s[1]);
	for(int i = M-3; i > 0; --i)
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

	Type *ta, *tb, *tc;
	ta = new Type[M+1];
	tb = new Type[M+1];
	ta[0] = Rod::t0[0];
	tb[0] = Rod::t0[0];
	for(int i = 1; i < M; ++i)
		ta[i] = Rod::t0[i];
	ta[M] = Rod::t0[M];
	tb[M] = Rod::t0[M];

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

	for(int j = 1; j < M; ++j)
	{
		tN[j] = ta[j];
	}

	delete [] K;
	delete [] ta;
	delete [] tb;
}

void TestCase4()
{
	Rod::Initialize(430);
	int NOI = Rod::M-1;

	float *s = new float[NOI];
	for(int i = 1; i < Rod::M; ++i)
	{
		s[i-1] = (float)i/(float)Rod::M;//(float)( (1 + tanh( (2*i/Rod::M - 1)*M_PI ))/2 );
	}

	Rod targetRod = Rod::CreateRod(s);
	Solve_HeatEquation(targetRod.s, targetRod.tN);
//	targetRod.Show();

	ADS::Initialize();
	ADV *S = new ADV[Rod::M-1];
	for(int i = 0; i < NOI; ++i)
	{
		S[i] = (float)0.5*s[i];
	}
	ADV *tN = new ADV[Rod::M+1];
	Solve_HeatEquation(S, tN);
	ADV sqe = SquaredError(NOI, targetRod.tN+1, tN+1);
	printf("SqE = %6.4f\n",sqe.v);
	float *J = new float[NOI];
	ADS::GetJacobianReverse(J, 1, NOI);
	ShowJacobian(J, 1, NOI, true);


	delete [] s;
	delete [] S;
	delete [] tN;
	targetRod.Destroy();
	ADS::Clear();
	Rod::Clear();
}