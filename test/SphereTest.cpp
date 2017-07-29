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

#include <terra/Sphere.hpp>
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
	TestContext() : sphere(6378137.0), tolerance(0.00001),
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
		{  1332415.577412, -4646808.068376,  4160833.916045 },
		{ -2509763.269724, -4645951.139178,  3578161.024896 },
		{ -1886194.019727, -5352272.552631, -2911590.264866 },
		{ -3951931.698399,  3351108.060135,  3719352.097806 },
		{  1845099.069584,  6106512.555871,   -76694.720684 }
	}
	{ }
	terra::Sphere<double> sphere;
	double tolerance;
	Coord<double>::type geod[6];
	Coord<double>::type ecef[6];
};
template<>
struct TestContext<float> {
	TestContext() : sphere(6378137.0f), tolerance(0.75f),
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
		{  1332415.577412f, -4646808.068376f,  4160833.916045f },
		{ -2509763.269724f, -4645951.139178f,  3578161.024896f },
		{ -1886194.019727f, -5352272.552631f, -2911590.264866f },
		{ -3951931.698399f,  3351108.060135f,  3719352.097806f },
		{  1845099.069584f,  6106512.555871f,   -76694.720684f }
	}
       	{}
	terra::Sphere<float> sphere;
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
testSphereSingleInplace(TestContext<T> const & ctx)
{
#define FUNC "testSphereSingleInplace: "
	for (auto i = 0u; i < sizeof ctx.geod/sizeof(typename Coord<T>::type); ++i) {
		typename Coord<T>::type coord = { ctx.geod[i][0], ctx.geod[i][1], ctx.geod[i][2] };
		terra::geodToECEF(&coord, ctx.sphere);
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

		terra::ecefToGeod(&coord, ctx.sphere);
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
testSphereSingle(TestContext<T> const &ctx)
{
#define FUNC "testSphereSingle: "
	for (auto i = 0u; i < sizeof ctx.geod/sizeof(typename Coord<T>::type); ++i) {
		typename Coord<T>::type ecef;
		terra::geodToECEF(&ecef, ctx.geod[i], ctx.sphere);
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
		terra::ecefToGeod(&geod, ecef, ctx.sphere);
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
testSphereSoA(TestContext<T> const &ctx)
{
#define FUNC "testSphereSoA: "

	constexpr auto const numCoords = sizeof ctx.geod/sizeof(typename Coord<T>::type);

	CoordSoA<T> geod, ecef;
	geod.x = new T[numCoords];
	geod.y = new T[numCoords];
	geod.z = new T[numCoords];
	ecef.x = new T[numCoords];
	ecef.y = new T[numCoords];
	ecef.z = new T[numCoords];

	createSoA(&geod, ctx.geod, numCoords);

	terra::geodToECEFSoA(&ecef, geod, numCoords, ctx.sphere);

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

	terra::ecefToGeodSoA(&geod, ecef, numCoords, ctx.sphere);

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
testSphereAoS(TestContext<T> const &ctx)
{
#define FUNC "testSphereSoA: "

	constexpr auto const numCoords = sizeof ctx.geod/sizeof(typename Coord<T>::type);

	typename Coord<T>::type* ecef = new typename Coord<T>::type[numCoords];
	typename Coord<T>::type* geod = new typename Coord<T>::type[numCoords];

	terra::geodToECEFAoS(&ecef, ctx.geod, numCoords, ctx.sphere);

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

	terra::ecefToGeodAoS(&geod, ecef, numCoords, ctx.sphere);

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
testSphere()
{
	TestContext<float> const ctxSP;
	TestContext<double> const ctxDP;

	testSphereSingleInplace(ctxSP);
	testSphereSingle(ctxSP);
	testSphereSoA(ctxSP);
	testSphereAoS(ctxSP);
	testSphereSingleInplace(ctxDP);
	testSphereSingle(ctxDP);
	testSphereSoA(ctxDP);
	testSphereAoS(ctxDP);
}
