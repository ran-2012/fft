
#pragma once

#include <algorithm>
#include <cstddef>
#include <complex>
#include <vector>
#include <cmath>

using namespace std;

constexpr double pi = 3.1415926;

auto raw_fft(const vector<double> &sig, const int type = 1)
{
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
		auto omega = [groupsize, type](double x)
		{
			++x;
			return exp(-complex<double>(0, type * 2 * x * pi / groupsize));
		};
		for (size_t n = 0; n < (size / groupsize); ++n)
		{
			for (size_t off = 0; off < r; ++off)
			{
				temp[n*groupsize + off] = ret[n*groupsize + off] + omega(off)*ret[n*groupsize + off + r];
				temp[n*groupsize + off + r] = ret[n*groupsize + off] - omega(off)*ret[n*groupsize + off + r];
			}
		}
		ret = move(temp);
	}
	return ret;
}

auto fft(const vector<double> &data)
{
	return raw_fft(data, 1);
}

auto afft(const vector<double> &data)
{
	return raw_fft(data, -1);
}