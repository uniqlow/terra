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

#ifndef terra_impl_SphereImpl_hpp
#define terra_impl_SphereImpl_hpp

namespace terra {

template<typename T, typename Coord>
inline
void
geodToECEF(
	Coord * const coord,
	Sphere<T> const sphere) noexcept
{
	assert(coord && "coord is nullptr");

	auto const lon = (*coord)[0];
	auto const lat = (*coord)[1];
	auto const alt = (*coord)[2];
	auto const n = sphere.radius + alt;
	auto const sin_lon = std::sin(lon);
	auto const cos_lon = std::cos(lon);
	auto const sin_lat = std::sin(lat);
	auto const cos_lat = std::cos(lat);
	(*coord)[0] = n*cos_lat*cos_lon;
	(*coord)[1] = n*cos_lat*sin_lon;
	(*coord)[2] = n*sin_lat;
}
 
template<typename T, typename Coord>
inline
void
geodToECEF(
	Coord * const TERRA_RESTRICT toECEF,
	Coord const & TERRA_RESTRICT fromGeodetic,
	Sphere<T> const sphere) noexcept
{
	assert(toECEF && "toECEF is nullptr");

	auto const lon = fromGeodetic[0];
	auto const lat = fromGeodetic[1];
	auto const alt = sphere.radius + fromGeodetic[2];
	auto const sin_lon = std::sin(lon);
	auto const cos_lon = std::cos(lon);
	auto const sin_lat = std::sin(lat);
	auto const cos_lat = std::cos(lat);
	(*toECEF)[0] = alt*cos_lat*cos_lon;
	(*toECEF)[1] = alt*cos_lat*sin_lon;
	(*toECEF)[2] = alt*sin_lat;
}

template<typename T, typename Coord>
inline
void
ecefToGeod(
	Coord * const coord,
	Sphere<T> const sphere) noexcept
{
	assert(coord && "coord is nullptr");

	auto const x = (*coord)[0];
	auto const y = (*coord)[1];
	auto const z = (*coord)[2];
	auto const r = sphere.radius;
	auto const p = std::sqrt(x*x + y*y);
	auto const lon = std::atan2(y, x);
	auto const lat = std::atan2(z, p);
	auto const coslat = std::cos(lat);
	auto const alt = (p/coslat) - r;
	(*coord)[0] = lon;
	(*coord)[1] = lat;
	(*coord)[2] = alt;
}

template<typename T, typename Coord>
inline
void
ecefToGeod(
	Coord * const TERRA_RESTRICT toGeodetic,
	Coord const & TERRA_RESTRICT fromECEF,
	Sphere<T> const sphere) noexcept
{
	assert(toGeodetic && "toGeodetic is nullptr");

	auto const r = sphere.radius;
	auto const x = fromECEF[0];
	auto const y = fromECEF[1];
	auto const z = fromECEF[2];
	auto const p = std::sqrt(x*x + y*y);
	auto const lon = std::atan2(y, x);
	auto const lat = std::atan2(z, p);
	auto const coslat = std::cos(lat);
	auto const alt = (p/coslat) - r;
	(*toGeodetic)[0] = lon;
	(*toGeodetic)[1] = lat;
	(*toGeodetic)[2] = alt;
}

template<typename T, typename Coord>
inline
void
geodToECEFSoA(
	Coord * const TERRA_RESTRICT toECEF,
	Coord const & TERRA_RESTRICT fromGeodetic,
	unsigned const numCoords,
	Sphere<T> const sphere) noexcept
{
	assert(toECEF && "toECEF is nullptr");

	auto const r = sphere.radius;

	for (auto i = 0u; i < numCoords; ++i) {
		auto const lon = fromGeodetic.x[i];
		auto const lat = fromGeodetic.y[i];
		auto const alt = r + fromGeodetic.z[i];
		auto const sin_lon = std::sin(lon);
		auto const cos_lon = std::cos(lon);
		auto const sin_lat = std::sin(lat);
		auto const cos_lat = std::cos(lat);
		toECEF->x[i] = alt*cos_lat*cos_lon;
		toECEF->y[i] = alt*cos_lat*sin_lon;
		toECEF->z[i] = alt*sin_lat;
	}	
}

template<typename T, typename Coord, typename Coord2>
inline
void
geodToECEFAoS(
	Coord * const TERRA_RESTRICT toECEF,
	Coord2 const & TERRA_RESTRICT fromGeodetic,
	unsigned const numCoords,
	Sphere<T> const sphere) noexcept
{
	assert(toECEF && "toECEF is nullptr");

	auto const r = sphere.radius;

	for (auto i = 0u; i < numCoords; ++i) {
		auto const lon = fromGeodetic[i][0];
		auto const lat = fromGeodetic[i][1];
		auto const alt = r + fromGeodetic[i][2];
		auto const sin_lon = std::sin(lon);
		auto const cos_lon = std::cos(lon);
		auto const sin_lat = std::sin(lat);
		auto const cos_lat = std::cos(lat);
		(*toECEF)[i][0] = alt*cos_lat*cos_lon;
		(*toECEF)[i][1] = alt*cos_lat*sin_lon;
		(*toECEF)[i][2] = alt*sin_lat;
	}	
}

template<typename T, typename Coord>
inline
void
ecefToGeodSoA(
	Coord * const TERRA_RESTRICT toGeodetic,
	Coord const & TERRA_RESTRICT fromECEF,
	unsigned const numCoords,
	Sphere<T> const sphere) noexcept
{
	assert(toGeodetic && "toGeodetic is nullptr");

	auto const r = sphere.radius;

	for (auto i = 0u; i < numCoords; ++i) {
		auto const x = fromECEF.x[i];
		auto const y = fromECEF.y[i];
		auto const z = fromECEF.z[i];
		auto const p = std::sqrt(x*x + y*y);
		auto const lon = std::atan2(y, x);
		auto const lat = std::atan2(z, p);
		auto const coslat = std::cos(lat);
		auto const alt = (p/coslat) - r;
		toGeodetic->x[i] = lon;
		toGeodetic->y[i] = lat;
		toGeodetic->z[i] = alt;
	}	
}

template<typename T, typename Coord, typename Coord2>
inline
void
ecefToGeodAoS(
	Coord * const TERRA_RESTRICT toGeodetic,
	Coord2 const & TERRA_RESTRICT fromECEF,
	unsigned const numCoords,
	Sphere<T> const sphere) noexcept
{
	assert(toGeodetic && "toGeodetic is nullptr");

	auto const r = sphere.radius;

	for (auto i = 0u; i < numCoords; ++i) {
		auto const x = fromECEF[i][0];
		auto const y = fromECEF[i][1];
		auto const z = fromECEF[i][2];
		auto const p = std::sqrt(x*x + y*y);
		auto const lon = std::atan2(y, x);
		auto const lat = std::atan2(z, p);
		auto const coslat = std::cos(lat);
		auto const alt = (p/coslat) - r;
		(*toGeodetic)[i][0] = lon;
		(*toGeodetic)[i][1] = lat;
		(*toGeodetic)[i][2] = alt;
	}	
}

} // !namespace terra

#endif // !terra_impl_SphereImpl_hpp
