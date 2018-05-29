
#include <fstream>
#include "fft.h"
#include <vector>
#include <complex>

using namespace std;
int main()
{
	ofstream out("data.txt");
	ofstream out2("raw.txt");
	//生成正弦波
	double pi = 3.1415926;
	auto generateSineWave = [pi](double f, double phi, double a)
	{
		return [=](double t)
		{
			return a * sin(2 * pi*f*t + phi);
		};
	};
	auto generateTriangleWave = [](double delta)
	{
		return [delta](double t)
		{
			if (abs(t) > delta)
				return (double)0.0;
			return 1 - abs(t) / delta;
		};
	};
	auto triangleWave = generateTriangleWave(0.2);
	//正弦波参数
	vector<vector<double>> sineWaveParameter =
	{
		{ 15,0,0.01 },
		{ 20,0,0.02 },
		{ 25,0,0.01 },
		{ 5,0,0.02 },
		{ 1,0,0.01 }
	};
	//信号
	auto signal = [&](double t)
	{
		double ret = 0;
		for (auto &i : sineWaveParameter)
			ret += generateSineWave(i[0], i[1], i[2])(t);
		return ret;
	};
	double delta = 0.01;
	int n = 1024;
	vector<double> raw;
	for (int i = 0; i < n; ++i)
	{
		raw.push_back(signal(i*delta));
	}
	auto freq = fft(raw);
	auto rawdata = afft(freq);
	for (auto i : freq)
	{
		out << abs(i) << endl;
	}
	for (int i = 0; i < n; ++i)
	{
		out2 << raw[i] << '\t' << rawdata[i].real() << endl;;
	}
	return 0;
}
