/*
 * Copyright (c) 2017 Jon Olsson <jlo@wintermute.net>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#include <terra/Ellipsoid.hpp>
#include <cstdio>
#include <cstdlib>

#define DEG2RAD(a) (a*(3.141592653/180.0f))

namespace {

template<typename T>
struct Coord {
	using type = T[3];
};

template<typename T>
struct CoordSoA {
	T* x;
	T* y;
	T* z;
};

template<typename T>
struct TestContext {
};
template<>
struct TestContext<double> {
	TestContext() : ellipsoid(6378137.0, 6356752.314245), tolerance(0.00001),
	geod{
		{ DEG2RAD(   0.000000), DEG2RAD(   0.000000),    0.0 },
		{ DEG2RAD( -74.000401), DEG2RAD(  40.719645),    5.0 },
		{ DEG2RAD(-118.378113), DEG2RAD(  34.122223),  500.0 },
		{ DEG2RAD(-109.412964), DEG2RAD( -27.160732),  100.0 },
		{ DEG2RAD( 139.703152), DEG2RAD(  35.671434),   50.0 },
		{ DEG2RAD(  73.187668), DEG2RAD(  -0.688815), 1500.0 }
	},
	ecef{
		{  6378137.000000,        0.000000,        0.000000 },
		{  1334317.624619, -4653441.470488,  4138879.637461 },
		{ -2512410.732611, -4650851.993119,  3557958.544077 },
		{ -1887510.983407, -5356009.574636, -2894118.577641 },
		{ -3956437.465399,  3354928.802309,  3698665.741603 },
		{  1845099.961938,  6106515.509193,   -76181.454642 }
	}
	{ }
	terra::Ellipsoid<double> ellipsoid;
	double tolerance;
	Coord<double>::type geod[6];
	Coord<double>::type ecef[6];
};
template<>
struct TestContext<float> {
	TestContext() : ellipsoid(6378137.0f, 6356752.314245f), tolerance(0.75f),
	geod{
		{ DEG2RAD(   0.000000f), DEG2RAD(   0.000000f),    0.0f },
		{ DEG2RAD( -74.000401f), DEG2RAD(  40.719645f),    5.0f },
		{ DEG2RAD(-118.378113f), DEG2RAD(  34.122223f),  500.0f },
		{ DEG2RAD(-109.412964f), DEG2RAD( -27.160732f),  100.0f },
		{ DEG2RAD( 139.703152f), DEG2RAD(  35.671434f),   50.0f },
		{ DEG2RAD(  73.187668f), DEG2RAD(  -0.688815f), 1500.0f }
	},
	ecef{
		{  6378137.000000f,        0.000000f,        0.000000f },
		{  1334317.624619f, -4653441.470488f,  4138879.637461f },
		{ -2512410.732611f, -4650851.993119f,  3557958.544077f },
		{ -1887510.983407f, -5356009.574636f, -2894118.577641f },
		{ -3956437.465399f,  3354928.802309f,  3698665.741603f },
		{  1845099.961938f,  6106515.509193f,   -76181.454642f }
	}
       	{}
	terra::Ellipsoid<float> ellipsoid;
	float tolerance;
	Coord<float>::type geod[6];
	Coord<float>::type ecef[6];
};

template<typename T>
struct Type {
};
template<>
struct Type<float> {
	static constexpr auto const str = "Float";
};
template<>
struct Type<double> {
	static constexpr auto const str = "Double";
};


template<typename T>
static
void
testEllipsoidSingleInplace(TestContext<T> const & ctx)
{
#define FUNC "testSllipsoidSingleInplace: "
	for (auto i = 0u; i < sizeof ctx.geod/sizeof(typename Coord<T>::type); ++i) {
		typename Coord<T>::type coord = { ctx.geod[i][0], ctx.geod[i][1], ctx.geod[i][2] };
		terra::geodToECEF(&coord, ctx.ellipsoid);
		if (std::abs(coord[0] - ctx.ecef[i][0]) > ctx.tolerance) {
			auto const diff = std::abs(coord[0] - ctx.ecef[i][0]);
			std::fprintf(stderr, FUNC "%s: FAIL: ECEF X coordinate failed: %f != %f, %f\n",
				     Type<T>::str, coord[0], ctx.ecef[i][0], diff);
			exit(-1);
		}
		if (std::abs(coord[1] - ctx.ecef[i][1]) > ctx.tolerance) {
			auto const diff = std::abs(coord[1] - ctx.ecef[i][1]);
			std::fprintf(stderr, FUNC "%s: FAIL: ECEF Y coordinate failed: %f != %f, %f\n",
				     Type<T>::str, coord[1], ctx.ecef[i][1], diff);
			exit(-1);
		}
		if (std::abs(coord[2] - ctx.ecef[i][2]) > ctx.tolerance) {
			auto const diff = std::abs(coord[2] - ctx.ecef[i][2]);
			std::fprintf(stderr, FUNC "%s: FAIL: ECEF Z coordinate failed: %f != %f, %f\n",
				     Type<T>::str, coord[2], ctx.ecef[i][2], diff);
			exit(-1);
		}

		terra::ecefToGeod(&coord, ctx.ellipsoid);
		if (std::abs(coord[0] - ctx.geod[i][0]) > ctx.tolerance) {
			auto const diff = std::abs(coord[0] - ctx.geod[i][0]);
			std::fprintf(stderr, FUNC "%s: FAIL: Geodetic X coordinate failed: %f != %f, %f\n",
				     Type<T>::str, coord[0], ctx.geod[i][0], diff);
			exit(-1);
		}
		if (std::abs(coord[1] - ctx.geod[i][1]) > ctx.tolerance) {
			auto const diff = std::abs(coord[1] - ctx.geod[i][1]);
			std::fprintf(stderr, FUNC "%s: FAIL: Geodetic Y coordinate failed: %f != %f, %f\n",
				     Type<T>::str, coord[1], ctx.geod[i][1], diff);
			exit(-1);
		}
		if (std::abs(coord[2] - ctx.geod[i][2]) > ctx.tolerance) {
			auto const diff = std::abs(coord[2] - ctx.geod[i][2]);
			std::fprintf(stderr, FUNC "%s: FAIL: Geodetic Y coordinate failed: %f != %f, %f\n",
				     Type<T>::str, coord[2], ctx.geod[i][2], diff);
			exit(-1);
		}
	}
	std::printf(FUNC "%s: SUCCESS\n", Type<T>::str);
#undef FUNC
}

template<typename T>
static
void
testEllipsoidSingle(TestContext<T> const &ctx)
{
#define FUNC "testEllipsoidSingle: "
	for (auto i = 0u; i < sizeof ctx.geod/sizeof(typename Coord<T>::type); ++i) {
		typename Coord<T>::type ecef;
		terra::geodToECEF(&ecef, ctx.geod[i], ctx.ellipsoid);
		if (std::abs(ecef[0] - ctx.ecef[i][0]) > ctx.tolerance) {
			auto const diff = std::abs(ecef[0] - ctx.ecef[i][0]);
			std::fprintf(stderr, FUNC "%s: FAIL: ECEF X coordinate failed: %f != %f, %f\n",
				     Type<T>::str, ecef[0], ctx.ecef[i][0], diff);
			exit(-1);
		}
		if (std::abs(ecef[1] - ctx.ecef[i][1]) > ctx.tolerance) {
			auto const diff = std::abs(ecef[1] - ctx.ecef[i][1]);
			std::fprintf(stderr, FUNC "%s: FAIL: ECEF Y coordinate failed: %f != %f, %f\n",
				     Type<T>::str, ecef[1], ctx.ecef[i][1], diff);
			exit(-1);
		}
		if (std::abs(ecef[2] - ctx.ecef[i][2]) > ctx.tolerance) {
			auto const diff = std::abs(ecef[2] - ctx.ecef[i][2]);
			std::fprintf(stderr, FUNC "%s: FAIL: ECEF Z coordinate failed: %f != %f, %f\n",
				     Type<T>::str, ecef[2], ctx.ecef[i][2], diff);
			exit(-1);
		}

		typename Coord<T>::type geod;
		terra::ecefToGeod(&geod, ecef, ctx.ellipsoid);
		if (std::abs(geod[0] - ctx.geod[i][0]) > ctx.tolerance) {
			auto const diff = std::abs(geod[0] - ctx.geod[i][0]);
			std::fprintf(stderr, FUNC "%s: FAIL: Geodetic X coordinate failed: %f != %f, %f\n",
				     Type<T>::str, geod[0], ctx.geod[i][0], diff);
			exit(-1);
		}
		if (std::abs(geod[1] - ctx.geod[i][1]) > ctx.tolerance) {
			auto const diff = std::abs(geod[1] - ctx.geod[i][1]);
			std::fprintf(stderr, FUNC "%s: FAIL: Geodetic Y coordinate failed: %f != %f, %f\n",
				     Type<T>::str, geod[1], ctx.geod[i][1], diff);
			exit(-1);
		}
		if (std::abs(geod[2] - ctx.geod[i][2]) > ctx.tolerance) {
			auto const diff = std::abs(geod[2] - ctx.geod[i][2]);
			std::fprintf(stderr, FUNC "%s: FAIL: Geodetic Y coordinate failed: %f != %f, %f\n",
				     Type<T>::str, geod[2], ctx.geod[i][2], diff);
			exit(-1);
		}
	}
	std::printf(FUNC "%s: SUCCESS\n", Type<T>::str);
#undef FUNC
}

template<typename T>
static
void
createSoA(CoordSoA<T>* soa, typename Coord<T>::type const* src, unsigned const numCoords)
{
	for (auto i = 0u; i < numCoords; ++i) {
		soa->x[i] = src[i][0];
		soa->y[i] = src[i][1];
		soa->z[i] = src[i][2];
	}
}

template<typename T>
static
void
testEllipsoidSoA(TestContext<T> const &ctx)
{
#define FUNC "testEllipsoidSoA: "

	constexpr auto const numCoords = sizeof ctx.geod/sizeof(typename Coord<T>::type);

	CoordSoA<T> geod, ecef;
	geod.x = new T[numCoords];
	geod.y = new T[numCoords];
	geod.z = new T[numCoords];
	ecef.x = new T[numCoords];
	ecef.y = new T[numCoords];
	ecef.z = new T[numCoords];

	createSoA(&geod, ctx.geod, numCoords);

	terra::geodToECEFSoA(&ecef, geod, numCoords, ctx.ellipsoid);

	for (auto i = 0u; i < numCoords; ++i) {
		if (std::abs(ecef.x[i] - ctx.ecef[i][0]) > ctx.tolerance) {
			auto const diff = std::abs(ecef.x[i] - ctx.ecef[i][0]);
			std::fprintf(stderr, FUNC "%s: FAIL: ECEF X coordinate failed: %f != %f, %f\n",
				     Type<T>::str, ecef.x[i], ctx.ecef[i][0], diff);
			exit(-1);
		}
		if (std::abs(ecef.y[i] - ctx.ecef[i][1]) > ctx.tolerance) {
			auto const diff = std::abs(ecef.y[i] - ctx.ecef[i][1]);
			std::fprintf(stderr, FUNC "%s: FAIL: ECEF Y coordinate failed: %f != %f, %f\n",
				     Type<T>::str, ecef.y[i], ctx.ecef[i][1], diff);
			exit(-1);
		}
		if (std::abs(ecef.z[i] - ctx.ecef[i][2]) > ctx.tolerance) {
			auto const diff = std::abs(ecef.z[i] - ctx.ecef[i][2]);
			std::fprintf(stderr, FUNC "%s: FAIL: ECEF Y coordinate failed: %f != %f, %f\n",
				     Type<T>::str, ecef.z[i], ctx.ecef[i][2], diff);
			exit(-1);
		}
	}

	terra::ecefToGeodSoA(&geod, ecef, numCoords, ctx.ellipsoid);

	for (auto i = 0u; i < numCoords; ++i) {
		if (std::abs(geod.x[i] - ctx.geod[i][0]) > ctx.tolerance) {
			auto const diff = std::abs(geod.x[i] - ctx.geod[i][0]);
			std::fprintf(stderr, FUNC "%s: FAIL: Geodetic X coordinate failed: %f != %f, %f\n",
				     Type<T>::str, geod.x[i], ctx.geod[i][0], diff);
			exit(-1);
		}
		if (std::abs(geod.y[i] - ctx.geod[i][1]) > ctx.tolerance) {
			auto const diff = std::abs(geod.y[i] - ctx.geod[i][1]);
			std::fprintf(stderr, FUNC "%s: FAIL: Geodetic Y coordinate failed: %f != %f, %f\n",
				     Type<T>::str, geod.y[i], ctx.geod[i][1], diff);
			exit(-1);
		}
		if (std::abs(geod.z[i] - ctx.geod[i][2]) > ctx.tolerance) {
			auto const diff = std::abs(geod.z[i] - ctx.geod[i][2]);
			std::fprintf(stderr, FUNC "%s: FAIL: Geodetic Y coordinate failed: %f != %f, %f\n",
				     Type<T>::str, geod.z[i], ctx.geod[i][2], diff);
			exit(-1);
		}
	}
	delete[] ecef.z;
	delete[] ecef.y;
	delete[] ecef.x;
	delete[] geod.z;
	delete[] geod.y;
	delete[] geod.x;

	std::printf(FUNC "%s: SUCCESS\n", Type<T>::str);
#undef FUNC
}

template<typename T>
static
void
testEllipsoidAoS(TestContext<T> const &ctx)
{
#define FUNC "testEllipsoidSoA: "

	constexpr auto const numCoords = sizeof ctx.geod/sizeof(typename Coord<T>::type);

	typename Coord<T>::type* ecef = new typename Coord<T>::type[numCoords];
	typename Coord<T>::type* geod = new typename Coord<T>::type[numCoords];

	terra::geodToECEFAoS(&ecef, ctx.geod, numCoords, ctx.ellipsoid);

	for (auto i = 0u; i < numCoords; ++i) {
		if (std::abs(ecef[i][0] - ctx.ecef[i][0]) > ctx.tolerance) {
			auto const diff = std::abs(ecef[i][0] - ctx.ecef[i][0]);
			std::fprintf(stderr, FUNC "%s: FAIL: ECEF X coordinate failed: %f != %f, %f\n",
				     Type<T>::str, ecef[i][0], ctx.ecef[i][0], diff);
			exit(-1);
		}
		if (std::abs(ecef[i][1] - ctx.ecef[i][1]) > ctx.tolerance) {
			auto const diff = std::abs(ecef[i][1] - ctx.ecef[i][1]);
			std::fprintf(stderr, FUNC "%s: FAIL: ECEF Y coordinate failed: %f != %f, %f\n",
				     Type<T>::str, ecef[i][1], ctx.ecef[i][1], diff);
			exit(-1);
		}
		if (std::abs(ecef[i][2] - ctx.ecef[i][2]) > ctx.tolerance) {
			auto const diff = std::abs(ecef[i][2] - ctx.ecef[i][2]);
			std::fprintf(stderr, FUNC "%s: FAIL: ECEF X coordinate failed: %f != %f, %f\n",
				     Type<T>::str, ecef[i][2], ctx.ecef[i][2], diff);
			exit(-1);
		}
	}

	terra::ecefToGeodAoS(&geod, ecef, numCoords, ctx.ellipsoid);

	for (auto i = 0u; i < numCoords; ++i) {
		if (std::abs(geod[i][0] - ctx.geod[i][0]) > ctx.tolerance) {
			auto const diff = std::abs(geod[i][0] - ctx.geod[i][0]);
			std::fprintf(stderr, FUNC "%s: FAIL: Geodetic X coordinate failed: %f != %f, %f\n",
				     Type<T>::str, geod[i][0], ctx.geod[i][0], diff);
			exit(-1);
		}
		if (std::abs(geod[i][1] - ctx.geod[i][1]) > ctx.tolerance) {
			auto const diff = std::abs(geod[i][1] - ctx.geod[i][1]);
			std::fprintf(stderr, FUNC "%s: FAIL: Geodetic Y coordinate failed: %f != %f, %f\n",
				     Type<T>::str, geod[i][1], ctx.geod[i][1], diff);
			exit(-1);
		}
		if (std::abs(geod[i][2] - ctx.geod[i][2]) > ctx.tolerance) {
			auto const diff = std::abs(geod[i][2] - ctx.geod[i][2]);
			std::fprintf(stderr, FUNC "%s: FAIL: Geodetic Y coordinate failed: %f != %f, %f\n",
				     Type<T>::str, geod[i][2], ctx.geod[i][2], diff);
			exit(-1);
		}
	}

	delete[] geod;
	delete[] ecef;

	std::printf(FUNC "%s: SUCCESS\n", Type<T>::str);
#undef FUNC
}

} // !namespace

void
testEllipsoid()
{
	TestContext<float> const ctxSP;
	TestContext<double> const ctxDP;

	testEllipsoidSingleInplace(ctxSP);
	testEllipsoidSingle(ctxSP);
	testEllipsoidSoA(ctxSP);
	testEllipsoidAoS(ctxSP);
	testEllipsoidSingleInplace(ctxDP);
	testEllipsoidSingle(ctxDP);
	testEllipsoidSoA(ctxDP);
	testEllipsoidAoS(ctxDP);
}
