
#pragma once

#include <algorithm>
#include <cstddef>
#include <complex>
#include <vector>
#include <cmath>

using namespace std;

template<typename datatype>
auto raw_fft(const vector<datatype> &sig, const int type = 1)
{
	constexpr double pi = 3.1415926;
	vector<complex<double>> ret;
	const auto size = sig.size();
	ret.resize(size);
	auto reserve=[](size_t a, size_t n)
	{
		size_t ret = 0;
		while (n>1)
		{
			ret <<= 1;
			ret += a & 0x1;
			a >>= 1;
			n >>= 1;
		}
		return ret;
	};
	
	for (int i = 0; i != size; ++i)
	{
		ret[reserve(i, size)] = sig[i];
	}
	for (size_t r = 1; r < size; r <<= 1)
	{
		vector<complex<double>> temp;
		temp.resize(size);
		auto groupsize = 2 * r;
		auto omega = [groupsize, type, pi](double x)
		{
			//++x;
			return exp(-complex<double>(0, type * 2 * x * pi / groupsize));
		};
		for (size_t k = 0; k < size; k+=groupsize)
		{
			for (size_t off = 0; off < r; ++off)
			{
				temp[k + off] 
					= ret[k + off] + omega(off)*ret[k + off + r];
				temp[k + off + r] 
					= ret[k + off] - omega(off)*ret[k + off + r];
			}
		}
		ret = move(temp);
	}
	return ret;
}

template<typename datatype>
auto fft(const vector<datatype> &data)
{
	return raw_fft(data, 1);
}

template<typename datatype>
auto afft(const vector<datatype> &data)
{
	auto ret = raw_fft(data, -1);
	for (auto &i : ret)
		i /= data.size();
	return ret;
}