#include "interpolation.hpp"
#include "basisfunctions.hpp"
#include "linalg.hpp"
#include <string>
#include <fstream>
#include <cmath>
#include <exception>

namespace cie
{
namespace splinekernel
{
    // Returns the control points for a b-spline curve with given degree that interpolates the given points on slide 183
    ControlPointsAndKnotVector interpolateWithBSplineCurve(const ControlPoints2D& interpolationPoints, size_t polynomialDegree)
    {

        size_t numberOfPoints = interpolationPoints[0].size();
        ControlPoints2D controlPoints;

        // Throw exception if number of x-values is not equal number of y-values
        if (interpolationPoints[0].size() != interpolationPoints[1].size())
        {
            throw interpolationPoints;
        }

        std::vector<double> t = centripetalParameterPositions(interpolationPoints);
        std::vector<double> knotVector = knotVectorUsingAveraging(t, polynomialDegree);

        // Initialize a matrix of shape functions
        linalg::Matrix N = linalg::Matrix(numberOfPoints, numberOfPoints, 0.0);

        for (size_t k = 0; k < numberOfPoints; ++k)
        {
            for (size_t i = 0; i < numberOfPoints; ++i)
            {
                N(k, i) = evaluateBSplineBasis(t[k], i, polynomialDegree, knotVector);

            }
        }

        // Solve linear system of equations
        controlPoints[0] = linalg::solve(N, interpolationPoints[0]);
        controlPoints[1] = linalg::solve(N, interpolationPoints[1]);


        return { {controlPoints},{knotVector} };
 
    }



    // compute the parameter positions ti for each Qi using the centripetal method on slide 177 
    std::vector<double> centripetalParameterPositions(const ControlPoints2D& interpolationPoints)
    {
        size_t numberOfPoints = interpolationPoints[0].size();
        std::vector<double> t(numberOfPoints, 0.0);

        t[0] = 0;
        t[numberOfPoints - 1] = 1;

        std::vector<double> d(numberOfPoints - 1, 0.0);
        double d_sum = 0;

        // get the sum of sqrt(d[i])(i = 0, 1, 2, ...., n-1)
        for (size_t i = 0; i < numberOfPoints - 1; ++i)
        {
            d[i] = sqrt(pow(interpolationPoints[0][i + 1] - interpolationPoints[0][i], 2) + pow(interpolationPoints[1][i + 1] - interpolationPoints[1][i], 2));

            d_sum = d_sum + sqrt(d[i]);
           
        }

        // get parameter postions t[i] (i = 1 ,2, 3, ..., n-1)
        for (size_t i = 0; i < numberOfPoints - 1; ++i)
        {
            t[i + 1] = t[i] + sqrt(d[i]) / d_sum;
        }

        return { t };
    }



    //compute the knot vector using the averaging technique on slide 180
    std::vector<double> knotVectorUsingAveraging(const std::vector<double>& parameterPositions,size_t polynomialDegree)
    {
        size_t numberOfPoints = parameterPositions.size();
        size_t numberOfKnots = numberOfPoints + polynomialDegree + 1;
        size_t numberOfInnerKnots = numberOfPoints - polynomialDegree - 1;
        if (polynomialDegree > numberOfPoints-1)
            throw std::runtime_error("Polynomial degree " + std::to_string(polynomialDegree) + " is too high for " + std::to_string(parameterPositions.size()) + " interpolation points.");
          
        std::vector<double> knotVector(numberOfKnots, 0.0);
        for (size_t i = 0; i < numberOfInnerKnots; ++i)
        {
            double sum = 0.0;
            for (size_t j = 1; j < polynomialDegree + 1; ++j)
            {
                sum += parameterPositions[i + j];
            }
            knotVector[i + polynomialDegree + 1] = sum / polynomialDegree;
        }
        for (size_t i = 0; i < polynomialDegree + 1; ++i)
        {
            knotVector[numberOfInnerKnots + polynomialDegree + 1 + i] = 1.0;
        }
        return knotVector;
       
    }



;

} // namespace splinekernel
} // namespace cie
