
// g++-4.8 -I../ThirdParty/ADOL-C/include -I../src/libcadet -std=c++0x -o compareSETFADadolC compareSETFADadolC.cpp

#include <iostream>
#include <type_traits>
#include <string>
#include <random>

#define ADOLC_TAPELESS
#define NUMBER_DIRECTIONS 4
#include <adolc/adouble.h>
ADOLC_TAPELESS_UNIQUE_INTERNALS

#define SFAD_DEFAULT_DIR 4
#include "setfad.hpp"
SFAD_GLOBAL_GRAD_SIZE

template <typename real_t> real_t sqr(const real_t& x) { return x * x; }

typedef sfad::FwdET<double, sfad::HeapStorage> FwdAD;

struct RandGen
{
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis;

    RandGen(double min, double max) : rd(), gen(rd()), dis(min, max) { }

    double operator()() { return dis(gen); }
};


inline bool isEqual(double a, double b) { return (std::abs(a-b) <= 1.0e-14) || (std::isnan(a) && std::isnan(b)); }


template <typename T>
void setADcombination(T& val, size_t combination)
{
	for (size_t i = 0; i < NUMBER_DIRECTIONS; ++i)
	{
		if (combination & (1 << i))
			val.setADValue(i, 1.0);
		else
			val.setADValue(i, 0.0);
	}
}


std::string combinationToString(size_t comb)
{
	std::string str(NUMBER_DIRECTIONS, '0');

	for (size_t i = 0; i < NUMBER_DIRECTIONS; ++i)
	{
		if (comb & (1 << i))
			str[i] = '1';
	}
	return str;
}


struct DataDependentTest { double data; };

template <typename T> struct testPlusL : public DataDependentTest { T operator()(const T& a) const { return a + data; } static const char* name() { return "testPlusL"; } };
template <typename T> struct testPlusR : public DataDependentTest { T operator()(const T& a) const { return data + a; } static const char* name() { return "testPlusR"; } };
template <typename T> struct testPlus { T operator()(const T& a, const T& b) const { return a + b; } static const char* name() { return "testPlus"; } };
template <typename T> struct testPlusE { T operator()(const T& a, const T& b) const { T c = a; c += b; return c; } static const char* name() { return "testPlusE"; } };

template <typename T> struct testMinusL : public DataDependentTest { T operator()(const T& a) const { return a - data; } static const char* name(){ return "testMinusL"; } };
template <typename T> struct testMinusR : public DataDependentTest { T operator()(const T& a) const { return data - a; } static const char* name(){ return "testMinusR"; } };
template <typename T> struct testMinus { T operator()(const T& a, const T& b) const { return a - b; } static const char* name() { return "testMinus"; } };
template <typename T> struct testMinusE { T operator()(const T& a, const T& b) const { T c = a; c -= b; return c; } static const char* name() { return "testMinusE"; } };

template <typename T> struct testMulL : public DataDependentTest { T operator()(const T& a) const { return a * data; } static const char* name() { return "testMulL"; } };
template <typename T> struct testMulR : public DataDependentTest { T operator()(const T& a) const { return data * a; } static const char* name() { return "testMulR"; } };
template <typename T> struct testMul { T operator()(const T& a, const T& b) const { return a * b; } static const char* name() { return "testMul"; } };
template <typename T> struct testMulE { T operator()(const T& a, const T& b) const { T c = a; c *= b; return c; } static const char* name() { return "testMulE"; } };

template <typename T> struct testDivL : public DataDependentTest { T operator()(const T& a) const { return a / data; } static const char* name() { return "testDivL"; } };
template <typename T> struct testDivR : public DataDependentTest { T operator()(const T& a) const { return data / a; } static const char* name() { return "testDivR"; } };
template <typename T> struct testDiv { T operator()(const T& a, const T& b) const { return a / b; } static const char* name() { return "testDiv"; } };
template <typename T> struct testDivE { T operator()(const T& a, const T& b) const { T c = a; c /= b; return c; } static const char* name() { return "testDivE"; } };

template <typename T> struct testUnaryPlus { T operator()(const T& a) const { return +a; } static const char* name() { return "testUnaryPlus"; } };
template <typename T> struct testUnaryMinus { T operator()(const T& a) const { return -a; } static const char* name() { return "testUnaryMinus"; } };


template<typename T> struct testExp { T operator()(const T &a) const { return exp(a); } static const char* name() { return "testExp"; } };
template<typename T> struct testLog { T operator()(const T &a) const { return log(a); } static const char* name() { return "testLog"; } };
template<typename T> struct testSqrt { T operator()(const T &a) const { return sqrt(a); } static const char* name() { return "testSqrt"; } };
template<typename T> struct testSqr { T operator()(const T &a) const { return sqr(a); } static const char* name() { return "testSqr"; } };

template<typename T> struct testSin { T operator()(const T &a) const { return sin(a); } static const char* name() { return "testSin"; } };
template<typename T> struct testCos { T operator()(const T &a) const { return cos(a); } static const char* name() { return "testCos"; } };
template<typename T> struct testTan { T operator()(const T &a) const { return tan(a); } static const char* name() { return "testTan"; } };
template<typename T> struct testAsin { T operator()(const T &a) const { return asin(a); } static const char* name() { return "testAsin"; } };
template<typename T> struct testAcos { T operator()(const T &a) const { return acos(a); } static const char* name() { return "testAcos"; } };
template<typename T> struct testAtan { T operator()(const T &a) const { return atan(a); } static const char* name() { return "testAtan"; } };

template<typename T> struct testPowB : public DataDependentTest { T operator()(const T &a) const { return pow(a, data); } static const char* name() { return "testPowB"; } };
template<typename T> struct testPowE : public DataDependentTest { T operator()(const T &a) const { return pow(data, a); } static const char* name() { return "testPowE"; } };
template<typename T> struct testPow { T operator()(const T &a, const T &b) const { return pow(a, b); } static const char* name() { return "testPow"; } };

template<typename T> struct testSinh { T operator()(const T &a) const { return sinh(a); } static const char* name() { return "testSinh"; } };
template<typename T> struct testCosh { T operator()(const T &a) const { return cosh(a); } static const char* name() { return "testCosh"; } };
template<typename T> struct testTanh { T operator()(const T &a) const { return tanh(a); } static const char* name() { return "testTanh"; } };

template<typename T> struct testFabs { T operator()(const T &a) const { return fabs(a); } static const char* name() { return "testFabs"; } };

template<typename T> struct testCeil { T operator()(const T &a) const { return ceil(a); } static const char* name() { return "testCeil"; } };
template<typename T> struct testFloor { T operator()(const T &a) const { return floor(a); } static const char* name() { return "testFloor"; } };

/*
template<typename T> struct testFmax { T operator()(const T &a, const FwdBase<T> &b) const { return fmax(a); } const char* name() { return "testFmax"; } };
template<typename T> struct testFmax { T operator()(T v, cTse<T> &a) const { return fmax(a); } const char* name() { return "testFmax"; } };
template<typename T> struct testFmax { T operator()(const T &a, T v) const { return fmax(a); } const char* name() { return "testFmax"; } };

template<typename T> struct testFmin { T operator()(const T &a, const FwdBase<T> &b) const { return fmin(a); } const char* name() { return "testFmin"; } };
template<typename T> struct testFmin { T operator()(T v, cTse<T> &a) const { return fmin(a); } const char* name() { return "testFmin"; } };
template<typename T> struct testFmin { T operator()(const T &a, T v) const { return fmin(a); } const char* name() { return "testFmin"; } };
*/


template<typename T, typename = void>
struct has_data : std::false_type { };

template<typename T>
struct has_data<T, decltype(std::declval<T>().data, void())> : std::true_type { };

template <class T, bool hasData = has_data<T>::value>
struct assignData { static void exec(T& f, double val)
{ 
	f.data = val;
} };

template <class T>
struct assignData<T, false> { static void exec(T&, double) { } };

template <class A, class B, bool hasData = has_data<A>::value && has_data<B>::value>
struct copyData { static void exec(const A& src, B& dest)
{
	dest.data = src.data;
} };

template <class A, class B>
struct copyData<A, B, false> { static void exec(A& src, B& dest) { } };


template <template<class T> class Func, class rnd_t>
void compareBinarySingle(rnd_t& rndGen)
{
	Func<adtl::adouble> f1;
	Func<FwdAD> f2;

	const double aVal = rndGen();
	const double bVal = rndGen();

	for (int j = 0; j < (1 << NUMBER_DIRECTIONS); ++j)
	{
		for (int k = 0; k < (1 << NUMBER_DIRECTIONS); ++k)
		{
			adtl::adouble a;
			a.setValue(aVal);
			setADcombination(a, j);

			adtl::adouble b;
			b.setValue(bVal);
			setADcombination(b, k);

			FwdAD c;
			c.setValue(aVal);
			setADcombination(c, j);

			FwdAD d;
			d.setValue(bVal);
			setADcombination(d, k);

			const adtl::adouble res1 = f1(a, b);
			const FwdAD res2 = f2(c, d);

			const std::string strCombA = combinationToString(j);
			const std::string strCombB = combinationToString(k);

			if (!isEqual(res1.getValue(), res2.getValue()))
				std::cout << "A = " << a.getValue() << " B = " << b.getValue() << " j = " << strCombA << " k = " << strCombB << ": Results differ " << res1.getValue() << " vs " << res2.getValue()
					<< " -> " << std::abs(res1.getValue() - res2.getValue()) << std::endl;

			for (int i = 0; i < NUMBER_DIRECTIONS; ++i)
			{
				if (!isEqual(res1.getADValue(i), res2.getADValue(i)))
					std::cout << "A = " << a.getValue() << " B = " << b.getValue() << " j = " << strCombA << " k = " << strCombB << ": Deriv " << i << " differs: " << res1.getADValue(i) << " vs " << res2.getADValue(i)
						<< " -> " << std::abs(res1.getADValue(i) - res2.getADValue(i)) << std::endl;
			}
		}	
	}
}


template <template<class T> class Func, class rnd_t>
void compareUnarySingle(rnd_t& rndGen)
{
	Func<adtl::adouble> f1;
	Func<FwdAD> f2;

	assignData< Func<adtl::adouble> >::exec(f1, rndGen());
	copyData< Func<adtl::adouble>, Func<FwdAD> >::exec(f1, f2);

	const double val = rndGen();

	for (size_t j = 0; j < (1 << NUMBER_DIRECTIONS); ++j)
	{
		adtl::adouble a;
		a.setValue(val);
		setADcombination(a, j);

		FwdAD c;
		c.setValue(val);
		setADcombination(c, j);

		const adtl::adouble res1 = f1(a);
		const FwdAD res2 = f2(c);

		const std::string strComb = combinationToString(j);

		if (!isEqual(res1.getValue(), res2.getValue()))
			std::cout << "Val = " << a.getValue() << " j = " << strComb << ": Results differ " << res1.getValue() << " vs " << res2.getValue()
				<< " -> " << std::abs(res1.getValue() - res2.getValue()) << std::endl;

		for (int i = 0; i < NUMBER_DIRECTIONS; ++i)
		{
			if (!isEqual(res1.getADValue(i), res2.getADValue(i)))
				std::cout << "Val = " << a.getValue() << " j = " << strComb << ": Deriv " << i << " differs: " << res1.getADValue(i) << " vs " << res2.getADValue(i)
					<< " -> " << std::abs(res1.getADValue(i) - res2.getADValue(i)) << std::endl;
		}
	}
}


template <template<class T> class Func, class rnd_t>
void compareBinary(rnd_t& rndGen)
{
	std::cout << "===== " << Func<adtl::adouble>::name() << std::endl;
	for (std::size_t i = 0; i < 20; ++i)
	{
		compareBinarySingle<Func, rnd_t>(rndGen);
	}
}


template <template<class T> class Func, class rnd_t>
void compareUnary(rnd_t& rndGen)
{
	std::cout << "===== " << Func<adtl::adouble>::name() << std::endl;
	for (std::size_t i = 0; i < 20; ++i)
	{
		compareUnarySingle<Func, rnd_t>(rndGen);
	}
}


int main(int argc, char** argv)
{
	RandGen dis(0.0, 10.0);

	compareBinary<testPlus>(dis);
	compareBinary<testPlusE>(dis);
	compareUnary<testPlusL>(dis);
	compareUnary<testPlusR>(dis);

	compareBinary<testMinus>(dis);
	compareBinary<testMinusE>(dis);
	compareUnary<testMinusL>(dis);
	compareUnary<testMinusR>(dis);

	compareBinary<testMul>(dis);
	compareBinary<testMulE>(dis);
	compareUnary<testMulL>(dis);
	compareUnary<testMulR>(dis);

	compareBinary<testDiv>(dis);
	compareBinary<testDivE>(dis);
	compareUnary<testDivL>(dis);
	compareUnary<testDivR>(dis);

	compareUnary<testUnaryPlus>(dis);
	compareUnary<testUnaryMinus>(dis);

	compareUnary<testExp>(dis);
	compareUnary<testLog>(dis);
	compareUnary<testSqrt>(dis);
	compareUnary<testSqr>(dis);
	compareUnary<testSin>(dis);
	compareUnary<testCos>(dis);
	compareUnary<testTan>(dis);
	compareUnary<testAsin>(dis);
	compareUnary<testAcos>(dis);
	compareUnary<testAtan>(dis);
	compareUnary<testPowB>(dis);
	compareUnary<testPowE>(dis);
	compareBinary<testPow>(dis);
	compareUnary<testSinh>(dis);
	compareUnary<testCosh>(dis);
	compareUnary<testTanh>(dis);
	compareUnary<testFabs>(dis);
	compareUnary<testCeil>(dis);
	compareUnary<testFloor>(dis);

	return 0;
}
