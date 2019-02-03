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

#ifndef terra_Sphere_hpp
#define terra_Sphere_hpp

#include <terra/Arch.hpp>
#include <cmath>
#include <cassert>

namespace terra {

/**
 * @brief Body approximated as a perfect sphere.
 * @tparam T floating-point type to be used (float or double).
 */
template<typename T>
struct Sphere {
	explicit Sphere(T const r) : radius(r) {}
	T radius;	/**< The radius of the sphere. */
};

/**
 * @brief Convert a geodetic coordinate to ECEF in place, using a reference sphere.
 * @tparam T floating-point type to be used (float or double).
 * @tparam Coord 3-tuple type that can be accessed via operator[].
 * @param coord Pointer to geodetic coordinate that will be overwritten by ECEF
 *	coordinate.
 *	Geodetic coordinate is indexed as: 0=longitude, 1=latitude, 2=altitude.
 *	ECEF coordinate is indexed as: 0=x, 1=y, 2=z.
 * @param sphere An instance of the reference sphere.
 */
template<typename T, typename Coord>
inline
void
geodToECEF(
	Coord * const coord,
	Sphere<T> const sphere) noexcept;

/**
 * @brief Convert a geodetic coordinate to an ECEF coordinate using a reference sphere.
 * @note: The two coordinates must not reference overlapping memory areas.
 * @tparam T floating-point type to be used (float or double).
 * @tparam Coord 3-tuple type that can be accessed via operator[].
 * @param toECEF Pointer to where the ECEF coordinate will be written.
 *	ECEF coordinate is indexed as: 0=x, 1=y, 2=z.
 * @param fromGeodetic The geodetic coordinate to be converted.
 *	Geodetic coordinate is indexed as: 0=longitude, 1=latitude, 2=altitude.
 * @param sphere An instance of the reference sphere.
 */
template<typename T, typename Coord>
inline
void
geodToECEF(
	Coord * const TERRA_RESTRICT toECEF,
	Coord const & TERRA_RESTRICT fromGeodetic,
	Sphere<T> const sphere) noexcept;

/**
 * @brief Convert an ECEF coordinate to geodetic in place, using a reference sphere.
 * @tparam T floating-point type to be used (float or double).
 * @tparam Coord 3-tuple type that can be accessed via operator[].
 * @param coord Pointer to ECEF coordinate that will be overwritten by geodetic
 *	coordinate.
 *	ECEF coordinate is indexed as: 0=x, 1=y, 2=z.
 *	Geodetic coordinate is indexed as: 0=longitude, 1=latitude, 2=altitude.
 * @param sphere An instance of the reference sphere.
 */
template<typename T, typename Coord>
inline
void
ecefToGeod(
	Coord * const coord,
	Sphere<T> const sphere) noexcept;

/**
 * @brief Convert an ECEF coordinate to a geodetic coordinate using a reference sphere.
 * @note: The two coordinates must not reference overlapping memory areas.
 * @tparam T floating-point type to be used (float or double).
 * @tparam Coord 3-tuple type that can be accessed via operator[].
 * @param toGeodetic Pointer to where the geodetic coordinate will be written.
 *	Geodetic coordinate is indexed as: 0=longitude, 1=latitude, 2=altitude.
 * @param fromECEF The ECEF coordinate to be converted.
 *	ECEF coordinate is indexed as: 0=x, 1=y, 2=z.
 * @param sphere An instance of the reference sphere.
 */
template<typename T, typename Coord>
inline
void
ecefToGeod(
	Coord * const TERRA_RESTRICT toGeodetic,
	Coord const & TERRA_RESTRICT fromECEF,
	Sphere<T> const sphere) noexcept;


/**
 * @brief Convert a series, in SoA form, of geodetic coordinates to ECEF coordinates using a reference sphere.
 * @note: The two coordinates must not reference overlapping memory areas.
 * @tparam T floating-point type to be used (float or double).
 * @tparam Coord a struct type with the arrays x, y, z for the coordinates.
 * @param toECEF Pointer to where the ECEF coordinates will be written.
 * @param fromGeodetic The geodetic coordinates to be converted.
 *	Geodetic coordinates are accessed as: x=longitude, y=latitude, z=altitude.
 * @param sphere An instance of the reference sphere.
 */
template<typename T, typename Coord>
inline
void
geodToECEFSoA(
	Coord * const TERRA_RESTRICT toECEF,
	Coord const & TERRA_RESTRICT fromGeodetic,
	unsigned const numCoords,
	Sphere<T> const sphere) noexcept;

/**
 * @brief Convert a series, in AoS form, of geodetic coordinates to ECEF coordinates using a reference sphere.
 * @note: The two coordinates must not reference overlapping memory areas.
 * @tparam T floating-point type to be used (float or double).
 * @tparam Coord 3-tuple type that can be accessed via operator[].
 * @tparam Coord2 3-tuple type that can be accessed via operator[].
 * @param toECEF Pointer to an array where the ECEF coordinates will be written.
 * @param fromGeodetic Pointer to an array of geodetic coordinates to be converted.
 *	Geodetic coordinates are indexed as: 0=longitude, 1=latitude, 2=altitude.
 * @param sphere An instance of the reference sphere.
 */
template<typename T, typename Coord, typename Coord2>
inline
void
geodToECEFAoS(
	Coord * const TERRA_RESTRICT toECEF,
	Coord2 const & TERRA_RESTRICT fromGeodetic,
	unsigned const numCoords,
	Sphere<T> const sphere) noexcept;

/**
 * @brief Convert a series, in SoA form, of ECEF coordinates to geodetic coordinates using a reference sphere.
 * @note: The two coordinates must not reference overlapping memory areas.
 * @tparam T floating-point type to be used (float or double).
 * @tparam Coord a struct type with the arrays x, y, z for the coordinates.
 * @param toGeodetic Pointer to where the geodetic coordinates will be written.
 *	Geodetic coordinates are accessed as: x=longitude, y=latitude, z=altitude.
 * @param fromECEF The ECEF coordinates to be converted.
 * @param sphere An instance of the reference sphere.
 */
template<typename T, typename Coord>
inline
void
ecefToGeodSoA(
	Coord * const TERRA_RESTRICT toGeodetic,
	Coord const & TERRA_RESTRICT fromECEF,
	unsigned const numCoords,
	Sphere<T> const sphere) noexcept;

/**
 * @brief Convert a series, in AoS form, of ECEF coordinates to geodetic coordinates using a reference sphere.
 * @note: The two coordinates must not reference overlapping memory areas.
 * @tparam T floating-point type to be used (float or double).
 * @tparam Coord 3-tuple type that can be accessed via operator[].
 * @tparam Coord2 3-tuple type that can be accessed via operator[].
 * @param toGeodetic Pointer to an array where the geodetic coordinates will be written.
 *	Geodetic coordinates are indexed as: 0=longitude, 1=latitude, 2=altitude.
 * @param fromECEF Pointer to an array of ECEF coordinates to be converted.
 * @param sphere An instance of the reference sphere.
 */
template<typename T, typename Coord, typename Coord2>
inline
void
ecefToGeodAoS(
	Coord * const TERRA_RESTRICT toGeodetic,
	Coord2 const * TERRA_RESTRICT fromECEF,
	unsigned const numCoords,
	Sphere<T> const sphere) noexcept;

} // !namespace terra

#include <terra/impl/SphereImpl.hpp>

#endif // !terra_Sphere_hpp
