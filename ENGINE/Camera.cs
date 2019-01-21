using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ENGINE
{
    public class Camera
    {
        public Vector3 Position;
        public Vector3 TargetPosition;
        public Vector3 Upvector;
        public DenseMatrix ViewMatrix;

        private Vector3 zAxis;
        private Vector3 xAxis;
        private Vector3 yAxis;
        private DenseMatrix CameraMatrix;

        public Camera(Vector3 position, Vector3 targetPosition)
        {
            Position = position;
            TargetPosition = targetPosition;
            Upvector = new Vector3(0, 0, 1);

            zAxis = Position - TargetPosition;
            zAxis.Normalize();

            xAxis = Upvector * zAxis;
            xAxis.Normalize();

            yAxis = zAxis * xAxis;
            yAxis.Normalize();

            CameraMatrix = DenseMatrix.OfArray(new double[,]
            {
                {xAxis.x, yAxis.x, zAxis.x, Position.x},
                {xAxis.y, yAxis.y, zAxis.y, Position.y},
                {xAxis.z, yAxis.z, zAxis.z, Position.z},
                { 0,0,0,1 }
            });

            ViewMatrix = DenseMatrix.OfMatrix(CameraMatrix.Inverse());
        }
    }
}
