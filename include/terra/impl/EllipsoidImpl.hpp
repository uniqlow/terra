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

#ifndef terra_impl_EllipsoidImpl_hpp
#define terra_impl_EllipsoidImpl_hpp

namespace terra {

template<typename T, typename Coord>
inline
void
geodToECEF(
	Coord * const coord,
	Ellipsoid<T> const ellipsoid) noexcept
{
	assert(coord && "coord is nullptr");

	auto const a = ellipsoid.major;
	auto const a2 = a*a;
	auto const b = ellipsoid.minor;
	auto const b2 = b*b;

	auto const lon = (*coord)[0];
	auto const lat = (*coord)[1];
	auto const alt = (*coord)[2];
	auto const sin_lon = std::sin(lon);
	auto const cos_lon = std::cos(lon);
	auto const sin_lat = std::sin(lat);
	auto const cos_lat = std::cos(lat);
	auto const Nphi = a2/(std::sqrt(a2*cos_lat*cos_lat + b2*sin_lat*sin_lat));
	auto const Nphi_alt_cos_lat = (Nphi + alt)*cos_lat;
	(*coord)[0] = Nphi_alt_cos_lat*cos_lon;
	(*coord)[1] = Nphi_alt_cos_lat*sin_lon;
	(*coord)[2] = ((b2/a2)*Nphi + alt)*sin_lat;
}

template<typename T, typename Coord>
inline
void
geodToECEF(
	Coord * const TERRA_RESTRICT toECEF,
	Coord const & TERRA_RESTRICT fromGeodetic,
	Ellipsoid<T> const ellipsoid) noexcept
{
	assert(toECEF && "toECEF is nullptr");

	auto const a = ellipsoid.major;
	auto const a2 = a*a;
	auto const b = ellipsoid.minor;
	auto const b2 = b*b;

	auto const lon = fromGeodetic[0];
	auto const lat = fromGeodetic[1];
	auto const alt = fromGeodetic[2];
	auto const sin_lon = std::sin(lon);
	auto const cos_lon = std::cos(lon);
	auto const sin_lat = std::sin(lat);
	auto const cos_lat = std::cos(lat);
	auto const Nphi = a2/(std::sqrt(a2*cos_lat*cos_lat + b2*sin_lat*sin_lat));
	auto const Nphi_alt_cos_lat = (Nphi + alt)*cos_lat;
	(*toECEF)[0] = Nphi_alt_cos_lat*cos_lon;
	(*toECEF)[1] = Nphi_alt_cos_lat*sin_lon;
	(*toECEF)[2] = ((b2/a2)*Nphi + alt)*sin_lat;
}

template<typename T, typename Coord>
inline
void
ecefToGeod(
	Coord * const coord,
	Ellipsoid<T> const ellipsoid) noexcept
{
	assert(coord && "coord is nullptr");

	auto const a = ellipsoid.major;
	auto const a2 = a*a;
	auto const b = ellipsoid.minor;
	auto const b2 = b*b;
	auto const e = std::sqrt((a2 - b2)/a2);
	auto const e2 = e*e;
	auto const ep = std::sqrt((a2 - b2)/b2);
	auto const ep2 = ep*ep;

	auto const x = (*coord)[0];
	auto const y = (*coord)[1];
	auto const z = (*coord)[2];

	auto const p = std::sqrt(x*x + y*y);
	auto const lon = std::atan2(y, x);
	auto const theta = std::atan2(z*a, p*b);
	auto const sin_theta = std::sin(theta);
	auto const cos_theta = std::cos(theta);
	auto const sin3_theta = sin_theta*sin_theta*sin_theta;
	auto const cos3_theta = cos_theta*cos_theta*cos_theta;
	auto const lat = std::atan2(z + ep2*b*sin3_theta, p - e2*a*cos3_theta);
	auto const cos_lat = std::cos(lat);
	auto const sin_lat = std::sin(lat);
	auto const N = a/(std::sqrt(1.0 - e2*sin_lat*sin_lat));
	auto const alt = (p/cos_lat) - N;
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
	Ellipsoid<T> const ellipsoid) noexcept
{
	assert(toGeodetic && "toGeodetic is nullptr");

	auto const a = ellipsoid.major;
	auto const a2 = a*a;
	auto const b = ellipsoid.minor;
	auto const b2 = b*b;
	auto const e = std::sqrt((a2 - b2)/a2);
	auto const e2 = e*e;
	auto const ep = std::sqrt((a2 - b2)/b2);
	auto const ep2 = ep*ep;

	auto const x = fromECEF[0];
	auto const y = fromECEF[1];
	auto const z = fromECEF[2];

	auto const p = std::sqrt(x*x + y*y);
	auto const lon = std::atan2(y, x);
	auto const theta = std::atan2(z*a, p*b);
	auto const sin_theta = std::sin(theta);
	auto const cos_theta = std::cos(theta);
	auto const sin3_theta = sin_theta*sin_theta*sin_theta;
	auto const cos3_theta = cos_theta*cos_theta*cos_theta;
	auto const lat = std::atan2(z + ep2*b*sin3_theta, p - e2*a*cos3_theta);
	auto const cos_lat = std::cos(lat);
	auto const sin_lat = std::sin(lat);
	auto const N = a/(std::sqrt(1.0 - e2*sin_lat*sin_lat));
	auto const alt = (p/cos_lat) - N;
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
	Ellipsoid<T> const ellipsoid) noexcept
{
	assert(toECEF && "toECEF is nullptr");

	auto const a = ellipsoid.major;
	auto const a2 = a*a;
	auto const b = ellipsoid.minor;
	auto const b2 = b*b;

	for (auto i = 0u; i < numCoords; ++i) {
		auto const lon = fromGeodetic.x[i];
		auto const lat = fromGeodetic.y[i];
		auto const alt = fromGeodetic.z[i];
		auto const sin_lon = std::sin(lon);
		auto const cos_lon = std::cos(lon);
		auto const sin_lat = std::sin(lat);
		auto const cos_lat = std::cos(lat);
		auto const Nphi = a2/(std::sqrt(a2*cos_lat*cos_lat + b2*sin_lat*sin_lat));
		auto const Nphi_alt_cos_lat = (Nphi + alt)*cos_lat;
		toECEF->x[i] = Nphi_alt_cos_lat*cos_lon;
		toECEF->y[i] = Nphi_alt_cos_lat*sin_lon;
		toECEF->z[i] = ((b2/a2)*Nphi + alt)*sin_lat;
	}
}

template<typename T, typename Coord, typename Coord2>
inline
void
geodToECEFAoS(
	Coord * const TERRA_RESTRICT toECEF,
	Coord2 const & TERRA_RESTRICT fromGeodetic,
	unsigned const numCoords,
	Ellipsoid<T> const ellipsoid) noexcept
{
	assert(toECEF && "toECEF is nullptr");

	auto const a = ellipsoid.major;
	auto const a2 = a*a;
	auto const b = ellipsoid.minor;
	auto const b2 = b*b;

	for (auto i = 0u; i < numCoords; ++i) {
		auto const lon = fromGeodetic[i][0];
		auto const lat = fromGeodetic[i][1];
		auto const alt = fromGeodetic[i][2];
		auto const sin_lon = std::sin(lon);
		auto const cos_lon = std::cos(lon);
		auto const sin_lat = std::sin(lat);
		auto const cos_lat = std::cos(lat);
		auto const Nphi = a2/(std::sqrt(a2*cos_lat*cos_lat + b2*sin_lat*sin_lat));
		auto const Nphi_alt_cos_lat = (Nphi + alt)*cos_lat;
		(*toECEF)[i][0] = Nphi_alt_cos_lat*cos_lon;
		(*toECEF)[i][1] = Nphi_alt_cos_lat*sin_lon;
		(*toECEF)[i][2] = ((b2/a2)*Nphi + alt)*sin_lat;
	}
}

template<typename T, typename Coord>
inline
void
ecefToGeodSoA(
	Coord * const TERRA_RESTRICT toGeodetic,
	Coord const & TERRA_RESTRICT fromECEF,
	unsigned const numCoords,
	Ellipsoid<T> const ellipsoid) noexcept
{
	assert(toGeodetic && "toGeodetic is nullptr");

	auto const a = ellipsoid.major;
	auto const a2 = a*a;
	auto const b = ellipsoid.minor;
	auto const b2 = b*b;
	auto const e = std::sqrt((a2 - b2)/a2);
	auto const e2 = e*e;
	auto const ep = std::sqrt((a2 - b2)/b2);
	auto const ep2 = ep*ep;

	for (auto i = 0u; i < numCoords; ++i) {
		auto const x = fromECEF.x[i];
		auto const y = fromECEF.y[i];
		auto const z = fromECEF.z[i];

		auto const p = std::sqrt(x*x + y*y);
		auto const lon = std::atan2(y, x);
		auto const theta = std::atan2(z*a, p*b);
		auto const sin_theta = std::sin(theta);
		auto const cos_theta = std::cos(theta);
		auto const sin3_theta = sin_theta*sin_theta*sin_theta;
		auto const cos3_theta = cos_theta*cos_theta*cos_theta;
		auto const lat = std::atan2(z + ep2*b*sin3_theta, p - e2*a*cos3_theta);
		auto const cos_lat = std::cos(lat);
		auto const sin_lat = std::sin(lat);
		auto const N = a/(std::sqrt(1.0 - e2*sin_lat*sin_lat));
		auto const alt = (p/cos_lat) - N;

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
	Ellipsoid<T> const ellipsoid) noexcept
{
	assert(toGeodetic && "toGeodetic is nullptr");

	auto const a = ellipsoid.major;
	auto const a2 = a*a;
	auto const b = ellipsoid.minor;
	auto const b2 = b*b;
	auto const e = std::sqrt((a2 - b2)/a2);
	auto const e2 = e*e;
	auto const ep = std::sqrt((a2 - b2)/b2);
	auto const ep2 = ep*ep;

	for (auto i = 0u; i < numCoords; ++i) {
		auto const x = fromECEF[i][0];
		auto const y = fromECEF[i][1];
		auto const z = fromECEF[i][2];

		auto const p = std::sqrt(x*x + y*y);
		auto const lon = std::atan2(y, x);
		auto const theta = std::atan2(z*a, p*b);
		auto const sin_theta = std::sin(theta);
		auto const cos_theta = std::cos(theta);
		auto const sin3_theta = sin_theta*sin_theta*sin_theta;
		auto const cos3_theta = cos_theta*cos_theta*cos_theta;
		auto const lat = std::atan2(z + ep2*b*sin3_theta, p - e2*a*cos3_theta);
		auto const cos_lat = std::cos(lat);
		auto const sin_lat = std::sin(lat);
		auto const N = a/(std::sqrt(1.0 - e2*sin_lat*sin_lat));
		auto const alt = (p/cos_lat) - N;
		(*toGeodetic)[i][0] = lon;
		(*toGeodetic)[i][1] = lat;
		(*toGeodetic)[i][2] = alt;
	}
}

} // !namespace terra

#endif // !terra_impl_EllipsoidImpl_hpp
