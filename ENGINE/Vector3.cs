using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System;

namespace ENGINE
{
    public class Vector3
    {
        public Matrix<double> vector;
        public double x => vector[0, 0];
        public double y => vector[1, 0];
        public double z => vector[2, 0];
        public double w => vector[3, 0];
        public double Length => Math.Sqrt(x * x + y * y + z * z);

        public Vector3(double X, double Y, double Z, double W = 1)
        {
            vector = DenseMatrix.OfArray(new double[,]
                {
                    {X},
                    {Y},
                    {Z},
                    {1}
                });
        }

        public static Vector3 operator +(Vector3 v1, Vector3 v2)
        {
            return new Vector3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
        }

        public static Vector3 operator -(Vector3 v1, Vector3 v2)
        {
            return new Vector3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
        }

        public static Vector3 operator /(Vector3 v, double a)
        {
            return new Vector3(v.x / a, v.y / a, v.z / a);
        }

        public static Vector3 operator *(Vector3 v, double a)
        {
            return new Vector3(v.x * a, v.y * a, v.z * a);
        }

        public static Vector3 operator *(Vector3 v1, Vector3 v2)
        {
            return new Vector3(v1.y*v2.z - v1.z*v2.y, v1.x*v2.z-v1.z*v2.x, v1.x*v2.y-v1.y*v2.x);
        }

        public static double Dot(Vector3 v1, Vector3 v2)
        {
            return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
        }

        public void Normalize()
        {
            var xn = x / Length;
            var yn = y / Length;
            var zn = z / Length;

            vector = DenseMatrix.OfArray(new double[,]
                {
                    {xn},
                    {yn},
                    {zn},
                    {w}
                });
        }
    }


}
