#pragma once

class Rod
{
public:
	static float L;
	static float T;
	static int M;
	static int N;
	static float lambda;
	static float *t0;

	float *s;
	float *tN;

	static void Initialize(int M);
	static void Clear();
	static Rod CreateRod(float *s =0, float *tN =0);
	void Destroy();
	void Show(bool isTerminal=true);
};



template <class Type>
void Solve_HeatEquation(Type *s, Type *tN);

