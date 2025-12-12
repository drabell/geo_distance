/*****************************************************************************************
Module          :  GeoDistance.cs | Class Lib | C#
Description     :  Methods to calculate the distance between two geo-points on Earth
Version         :  8.1.001
*****************************************************************************************
Author          :  Alexander Bell
Copyright       :  2011-2025 Alexander Bell
*****************************************************************************************
DISCLAIMER      :  This Module is provided on AS IS basis without any warranty.
                :  The user assumes the entire risk as to the accuracy and the use
                :  of this module. In no event shall the author be liable for any damages
                :  arising out of the use of or inability to use this module.
TERMS OF USE    :  This module is copyrighted. Please keep the Copyright notice intact.
*****************************************************************************************/
using System;

namespace GIS
{
    /// <summary>
    /// Class GeoDistance contains four static methods to calculate the
    /// great-circle (orthodromic) distance between two geo-points on Earth
    /// specified by coordinates in decimal format (Latitude, Longitude), e.g.
    /// John F. Kennedy International Airport (JFK): {40.641766,-73.780968},
    /// Los Angeles International Airport (LAX): {33.942791,-118.410042}
    /// Sample output:
    /// ===============================================================================
    /// Great-circle (orthodromic) distance between two geo-points:
    /// JFK {40.641766,-73.780968} to LHR {51.470020,-0.454295}
    /// km ----------------------------------------------------------------------------
    /// Haversine					: 5540.175419079548 (high accuracy)
    /// Spherical Law of Cosines	: 5540.175419079548 (high accuracy)
    /// Inverse Vincenty			: 5555.065686009474 (highest accuracy)
    /// Spherical Earth Projection	: 5784.908563389233 (lower accuracy)
    /// Expected value              :~5554.5 km
    /// miles -------------------------------------------------------------------------
    /// Haversine					: 3442.5054053574304 (high accuracy)
    /// Spherical Law of Cosines	: 3442.5054053574304 (high accuracy)
    /// Inverse Vincenty			: 3451.7577882724104 (highest accuracy)
    /// Spherical Earth Projection	: 3594.5755310171303 (lower accuracy)
    /// Expected value              :~3451.4 miles
    /// ===============================================================================
    /// </summary>
    public static class GeoDistance
    {
        // SI: km, US: miles
        public enum UnitSystem { SI = 0, US = 1 }

        #region private: const
        // Earth mean radius, km
        private const double _meanRadius = 6371.009;
        // Conversion factor: mile to km
        private const double _mi2km = 1.609344;
        // Conversion factor: degree to radian
        private const double _toRad = Math.PI / 180.0;
        #endregion

        #region Haversine algorithm ****************************************************
        /// <summary>
        /// Haversine algorithm implemented enables high-accuracy geodesic calculation 
        /// of the great-circle (a.k.a. orthodromic) distance (km/miles) between two 
        /// geographic points on the Earth's surface.
        /// </summary>
        /// <param name="Lat1">double: 1st point Latitude</param>
        /// <param name="Lon1">double: 1st point Longitude</param>
        /// <param name="Lat2">double: 2nd point Latitude</param>
        /// <param name="Lon2">double: 2nd point Longitude</param>
        /// <returns>double: distance, km/miles</returns>
        public static double Haversine(double Lat1, double Lon1,
                                       double Lat2,  double Lon2,
                                       UnitSystem UnitSys){
            try {
                double φ1 = Lat1 * _toRad; // Lat1 in radians;
                double φ2 = Lat2 * _toRad; // Lat2 in radians;

                double _a = Math.Sin((φ2 - φ1) / 2);
                _a *= _a; // calculate square

                double _b = Math.Sin(((Lon2 - Lon1)/2) *_toRad);
                _b *= _b * Math.Cos(φ1) * Math.Cos(φ2);

                // central angle, a.k.a. arc segment angular distance
                double _ca = 2 * Math.Asin(Math.Sqrt(_a + _b));

                // orthodromic distance on Earth between 2 points, km or miles
                return _ca * (UnitSys == UnitSystem.SI ? 1 : 1 / _mi2km) * _meanRadius;
            }
            catch { return -1; } //indicates error
        }
        #endregion

        #region Spherical Law of Cosines algorithm *************************************
        /// <summary>
        /// Spherical Law of Cosines (SLC) algorithm implemented in this method enables
        /// high-accuracy geodesic calculation of the great-circle (a.k.a. orthodromic) 
        /// distance (km/miles) between two geographic points on the Earth's surface.
        /// Note: results are very close to the Haversine formula, which is generally 
        /// preferred for numerical stability with small distances calculation.
        /// </summary>
        /// <param name="Lat1">double: 1st point Latitude</param>
        /// <param name="Lon1">double: 1st point Longitude</param>
        /// <param name="Lat2">double: 2nd point Latitude</param>
        /// <param name="Lon2">double: 2nd point Longitude</param>
        /// <returns>double: distance, km/miles</returns>
        public static double SLC(double Lat1, double Lon1,
                                 double Lat2, double Lon2,
                                 UnitSystem UnitSys){
            try {

                double φ1 = Lat1 * _toRad; // Lat1;
                double φ2 = Lat2 * _toRad; // Lat2;
                double Δλ = (Lon1 - Lon2) * _toRad;

                // central angle, aka arc segment angular distance
                double _ca = Math.Acos(Math.Sin(φ1) * Math.Sin(φ2) +
                        Math.Cos(φ1) * Math.Cos(φ2) * Math.Cos(Δλ));

                // orthodromic distance on Earth between 2 points, km or miles
                return (UnitSys == UnitSystem.SI ? 1 : 1 / _mi2km) * _ca * _meanRadius;
            }
            catch { return -1; } //indicates error
        }
        #endregion

        //region Vincenty algorithm (high accuracy) *************************************
        /// <summary>
        /// Inverse Vincenty (ellipsoid) algorithm implemented in this method enables
        /// the very high-accuracy geodesic calculation of the great-circle  (orthodromic)
        /// distance (km/miles) between two geographic points on the Earth's surface.
        /// Notes -----------------------------------------------------------------------
        /// Inverse Vincenty (ellipsoid) algorithm provides the highest accuracy among
        /// the common spherical/ellipsoidal computational methods, but it is not a 
        /// closed-form. This inverse solution (distance and bearings between two points)
        /// is an efficient iterative algorithm with nested expressions well-suited for
        /// the software implementation. Regarding its accuracy and robustness:
        /// - Convergence:
        /// The inverse method can fail near antipodal points.
        /// Use a max-iteration guard and a small epsilon; if it fails, fall back
        /// to a more robust geodesic algorithm.
        /// - Precision:
        /// Double precision is sufficient; avoid premature rounding of inputs.
        /// Keep lat/lon in radians for the loop.
        /// - Model choice:
        /// WGS84 is standard. For different datum (e.g., GRS80), set 𝑎/𝑓 accordingly.
        /// - Outputs:
        /// Besides distance, this method can return initial/final bearings.
        /// - AI vibe coding:
        /// This Inverse Vincenty geodesic algorithm was implemented in AI-assisted
        /// pair programming (vibe coding) interactive session with AI Copilot.
        /// -----------------------------------------------------------------------------
        /// </summary>
        /// <returns>double: orthodromic distance, km/miles</returns>
        public static double Vincenty(double lat1, double lon1,
                                      double lat2, double lon2,
                                      UnitSystem UnitSys)
        {
            // WGS84 constants
            double a = 6378137.0; // Earth equatorial radius, m
            double f = 1.0 / 298.257223563;
            double b = a * (1.0 - f);
            try
            {
                // Convert to radians
                double φ1 = (lat1) * _toRad, φ2 = (lat2) * _toRad;
                double Δλ = (lon2 - lon1) * _toRad;

                // Reduced latitudes
                double U1 = Math.Atan((1 - f) * Math.Tan(φ1));
                double U2 = Math.Atan((1 - f) * Math.Tan(φ2));

                double sinU1 = Math.Sin(U1), cosU1 = Math.Cos(U1);
                double sinU2 = Math.Sin(U2), cosU2 = Math.Cos(U2);

                double λ = Δλ;
                double λPrev;
                double iterLimit = 100;
                double ε = 1e-12;

                double sinσ, cosσ, σ, sinα, cos2α, cos2σm;
                double u2, A, B, Δσ;
                do
                {
                    double sinλ = Math.Sin(λ), cosλ = Math.Cos(λ);
                    double term1 = cosU2 * sinλ;
                    double term2 = cosU1 * sinU2 - sinU1 * cosU2 * cosλ;

                    sinσ = Math.Sqrt(term1 * term1 + term2 * term2);
                    if (sinσ == 0.0) return 0.0; // coincident points

                    cosσ = sinU1 * sinU2 + cosU1 * cosU2 * cosλ;
                    σ = Math.Atan2(sinσ, cosσ);

                    sinα = (cosU1 * cosU2 * sinλ) / sinσ;
                    double sin2α = sinα * sinα;
                    cos2α = 1.0 - sin2α;

                    if (cos2α != 0.0) cos2σm = cosσ - (2.0 * sinU1 * sinU2) / cos2α;
                    else cos2σm = 0.0; // equatorial line

                    u2 = (cos2α * (a * a - b * b)) / (b * b);

                    A = 1.0 + (u2 / 16384.0) * 
                        (4096.0 + u2 * (-768.0 + u2 * (320.0 - 175.0 * u2)));
                    B = (u2 / 1024.0) * 
                        (256.0 + u2 * (-128.0 + u2 * (74.0 - 47.0 * u2)));

                    double cos2σm2 = cos2σm * cos2σm;
                    Δσ = B * sinσ * (cos2σm + (B / 4.0) * (cosσ * (-1.0 + 2.0 * cos2σm2)
                            - (B / 6.0) * cos2σm * (-3.0 + 4.0 * sinσ * sinσ) * 
                            (-3.0 + 4.0 * cos2σm2)));

                    double C = (f / 16.0) * cos2α * (4.0 + f * (4.0 - 3.0 * cos2α));

                    λPrev = λ;
                    λ = Δλ + (1.0 - C) * f * sinα * (σ + C * sinσ * 
                        (cos2σm + C * cosσ * (-1.0 + 2.0 * cos2σm2)));

                    if (Math.Abs(λ - λPrev) < ε) break;
                } while (--iterLimit > 0);

                // If not converged, try to fall back to a robust algorithm here
                if (iterLimit == 0) throw new ArithmeticException("No Convergence");

                double s = b * A * (σ - Δσ);

                // Optional: initial/final bearings
                // double α1 = Math.atan2(cosU2 * Math.sin(λ),
                // cosU1 * sinU2 - sinU1 * cosU2 * Math.cos(λ));
                // double α2 = Math.atan2(cosU1 * Math.sin(λ),
                // -sinU1 * cosU2 + cosU1 * sinU2 * Math.cos(λ));

                // orthodromic distance on Earth between 2 points, km or miles
                return s * (UnitSys == UnitSystem.SI ? 1 : 1 / _mi2km) / 1000;
            }
            catch (Exception e) { return -1; } //indicates error
        }
        //endregion

        #region Spherical Earth Projection algorithm ***********************************
        /// <summary>
        /// Spherical Earth Projection (SEP) to a plane formula
        /// implemented in this method enables the calculation
        /// of a great-circle (orthodromic) distance(km/miles) between two
        /// geographic points on the Earth using Pythagorean Theorem:
        /// Central Angle: a = Sqrt((φ2 - φ1)^2 + (Cos((φ1 + φ2)/2) * (Lon2 - Lon1))^2)
        /// ---------------------------------------------------------------------------
        /// Note: this is a relatively low accuracy computation approach
        /// suitable for small distances (e.g., within a city or small region);
        /// it is shown mostly for a didactic purpose. For higher accuracy over
        /// longer distances, use either Haversine, or Spherical Law of Cosines,
        /// or Inverse Vincenty methods (the latter provides the highest accuracy).
        /// </summary>
        /// <param name="Lat1">double: 1st point Latitude</param>
        /// <param name="Lon1">double: 1st point Longitude</param>
        /// <param name="Lat2">double: 2nd point Latitude</param>
        /// <param name="Lon2">double: 2nd point Longitude</param>
        /// <returns>double: distance, km/miles</returns>
        public static double SEP(double Lat1, double Lon1,
                                 double Lat2, double Lon2,
                                 UnitSystem UnitSys){
            try {
                double φ1 = Lat1 * _toRad;
                double φ2 = Lat2 * _toRad;
                double Δλ = (Lat2 - Lat1) * _toRad;
                
                double _a = (Lon2 - Lon1) * Math.Cos((φ1 + φ2) / 2) * _toRad;

                // central angle, a.k.a. arc segment angular distance
                double _ca = Math.Sqrt(_a * _a + Δλ * Δλ);

                // orthodromic distance on Earth between 2 points, km or miles
                return _ca * (UnitSys == UnitSystem.SI ? 1 : 1 / _mi2km) * _meanRadius;
            }
            catch { return -1; } //indicates error
        }
        #endregion
    }
}
