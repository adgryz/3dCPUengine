using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections.Generic;
using System.Numerics;
using System.Runtime.InteropServices.WindowsRuntime;
using System.Threading.Tasks;
using Windows.Foundation;
using Windows.Storage;
using Windows.UI;
using Windows.UI.Xaml.Media.Imaging;
using static ENGINE.EnumModes;

namespace ENGINE
{


    class Device
    {
        private byte[] backBuffer;
        private readonly double[] zBuffer;
        private object[] lockBuffer;

        private WriteableBitmap bmp;
        private int screenWidth;
        private int ScreenHeight;

        double alpha = 30;
        double kd = 0.5;
        double ks = 1;
        Vector3 LP = new Vector3(0, 10 , 10);
        public ShadingMode shadingMode;
        public ReflectionMode reflectionMode;

        public Device(WriteableBitmap bmp, ShadingMode sMode, ReflectionMode rMode)
        {
            this.bmp = bmp;
            screenWidth = bmp.PixelWidth;
            ScreenHeight = bmp.PixelHeight;
            backBuffer = new byte[screenWidth * ScreenHeight * 4];
            zBuffer = new double[screenWidth * ScreenHeight];
            lockBuffer = new object[screenWidth * ScreenHeight];
            shadingMode = sMode;
            reflectionMode = rMode;
            for (var i = 0; i < lockBuffer.Length; i++)
            {
                lockBuffer[i] = new object();
            }
        }

        public void Clear(byte r, byte g, byte b, byte a)
        {
            for (var index = 0; index < backBuffer.Length; index += 4)
            {
                backBuffer[index] = b;
                backBuffer[index + 1] = g;
                backBuffer[index + 2] = r;
                backBuffer[index + 3] = a;
            }

            for (var index = 0; index < zBuffer.Length; index++)
            {
                zBuffer[index] = double.MaxValue;
            }
        }

        public void DrawScene()
        {
            using (var stream = bmp.PixelBuffer.AsStream())
            {
                stream.Write(backBuffer, 0, backBuffer.Length);
            }
            bmp.Invalidate();
        }

        public void PutPixel(int x, int y,double z, Color color)
        {
            var zIndex = x + y * screenWidth;
            var index = zIndex * 4;

            lock (lockBuffer[zIndex])
            {
                if (zBuffer[zIndex] < z)
                {
                    return;
                }

                zBuffer[zIndex] = z;

                backBuffer[index] = color.B;
                backBuffer[index + 1] = color.G;
                backBuffer[index + 2] = color.R;
                backBuffer[index + 3] = color.A;
            }
        }

        public void DrawPoint(Vector3 point, Color color)
        {
            if (point.x >= 0 && point.y >= 0 && point.x < screenWidth && point.y < ScreenHeight)
            {
                PutPixel((int)point.x, (int)point.y, point.z, color);
            }
        }

        public Vertex ProjectVertex(Vertex vertex, DenseMatrix transformMatrix, DenseMatrix worldMatrix)
        {
            var point2d = transformMatrix * vertex.Position.vector;
            var point3dWorld = worldMatrix * vertex.Position.vector;
            var normal3dWorld = worldMatrix * vertex.Normal.vector;

            var x = (point2d[0, 0] / point2d[3, 0] + 5)*50;
            var y = (point2d[1, 0] / point2d[3, 0] + 5)*50;
            var newV =  new Vertex
            {
                Position = new Vector3(x,y,point2d[2,0]),
                WorldPosition = new Vector3(point3dWorld[0,0] ,point3dWorld[1,0], point3dWorld[2,0]),
                Normal = new Vector3(normal3dWorld[0,0], normal3dWorld[1, 0], normal3dWorld[2, 0])
            };

            return newV;
        }

        public void Render(double Time, Camera camera, params Mesh[] meshes)
        {
            var ViewMatrix = camera.ViewMatrix;

            var ProjectionMatrix = DenseMatrix.OfArray(new double[,]
                {
                    {2.414,0,0,0},
                    {0,2.414,0,0},
                    {0,0,-1.02,-2.02},
                    {0,0,-1,0},
                });

            foreach (Mesh mesh in meshes)
            {
                var transformMatrix = ProjectionMatrix * ViewMatrix * mesh.ModelMatrix;

                Parallel.For(0, mesh.Faces.Length, faceIndex =>
                {
                    var face = mesh.Faces[faceIndex];
                    var vertexA = mesh.Vertices[face.A];
                    var vertexB = mesh.Vertices[face.B];
                    var vertexC = mesh.Vertices[face.C];

                    var projVerA = ProjectVertex(vertexA, transformMatrix, mesh.ModelMatrix);
                    var projVerB = ProjectVertex(vertexB, transformMatrix, mesh.ModelMatrix);
                    var projVerC = ProjectVertex(vertexC, transformMatrix, mesh.ModelMatrix);

                    var color = 255;
                    var isBlin = reflectionMode == ReflectionMode.Blin;

                    switch(shadingMode)
                    {
                        case ShadingMode.Flat:
                            DrawTriangle_FS(projVerA, projVerB, projVerC, Color.FromArgb(255, Limit(color), Limit(color), Limit(color)), Time, camera.Position, isBlin);
                            break;
                        case ShadingMode.Goraud:
                            DrawTriangle_GS(projVerA, projVerB, projVerC, Color.FromArgb(255, Limit(color), Limit(color), Limit(color)), Time, camera.Position, isBlin);
                            break;
                        case ShadingMode.Phong:
                            DrawTriangle_PS(projVerA, projVerB, projVerC, Color.FromArgb(255, Limit(color), Limit(color), Limit(color)), Time, camera.Position, isBlin);
                            break;
                    }

                    faceIndex++;
                });

            }
        }

        public async Task<Mesh[]> LoadJSONFileAsync(string fileName)
        {
            var meshes = new List<Mesh>();

            StorageFile file = await StorageFile.GetFileFromApplicationUriAsync(new Uri(fileName));
            var data = await Windows.Storage.FileIO.ReadTextAsync(file);
            dynamic jsonObject = Newtonsoft.Json.JsonConvert.DeserializeObject(data);

            for (var meshIndex = 0; meshIndex < jsonObject.meshes.Count; meshIndex++)
            {
                var verticesArray = jsonObject.meshes[meshIndex].vertices;
                var indicesArray = jsonObject.meshes[meshIndex].indices;
                var uvCount = jsonObject.meshes[meshIndex].uvCount.Value;
                var verticesStep = 1;

                switch ((int)uvCount)
                {
                    case 0:
                        verticesStep = 6;
                        break;
                    case 1:
                        verticesStep = 8;
                        break;
                    case 2:
                        verticesStep = 10;
                        break;
                }

                var verticesCount = verticesArray.Count / verticesStep;
                var facesCount = indicesArray.Count / 3;
                var mesh = new Mesh(jsonObject.meshes[meshIndex].name.Value, verticesCount, facesCount);

                for (var index = 0; index < verticesCount; index++)
                {
                    var x = (double)verticesArray[index * verticesStep].Value;
                    var y = (double)verticesArray[index * verticesStep + 1].Value;
                    var z = (double)verticesArray[index * verticesStep + 2].Value;

                    var nx = (double)verticesArray[index * verticesStep + 3].Value;
                    var ny = (double)verticesArray[index * verticesStep + 4].Value;
                    var nz = (double)verticesArray[index * verticesStep + 5].Value;

                    mesh.Vertices[index] = new Vertex
                    {
                        Position = new Vector3(x, y, z),
                        Normal = new Vector3(nx, ny, nz)
                    };
                }

                for (var index = 0; index < facesCount; index++)
                {
                    var a = (int)indicesArray[index * 3].Value;
                    var b = (int)indicesArray[index * 3 + 1].Value;
                    var c = (int)indicesArray[index * 3 + 2].Value;
                    mesh.Faces[index] = new Face { A = a, B = b, C = c };
                }

                var position = jsonObject.meshes[meshIndex].position;
                mesh.ModelMatrix = DenseMatrix.OfArray(new double[,]
                {
                        {position[0].Value,0,0,0},
                        {0,position[1].Value,0,0},
                        {0,0,position[2].Value,0},
                        {0,0,0,1},
                });
                meshes.Add(mesh);
            }
            return meshes.ToArray();
        }

        #region TriangleDrawing
        double Clamp(double value, double min = 0, double max = 1)
        {
            return Math.Max(min, Math.Min(value, max));
        }

        double Interpolate(double min, double max, double gradient)
        {
            return min + (max - min) * Clamp(gradient);
        }

        Vector3 InterpolateVector(Vector3 min, Vector3 max, double gradient)
        {
            return min + (max - min) * Clamp(gradient);
        }

        private void Swap(ref Vertex x, ref Vertex y)
        {
            var temp = y;
            y = x;
            x = temp;
        }

        private double OverZero(double value)
        {
            return Math.Max(0, value);
        }

        private byte Limit(double value)
        {
            return value > 255 ? (byte)255 : (byte)value;
        }
        #endregion

        double Compute_NoL(Vector3 point, Vector3 normal, Vector3 lightPosition)
        {
            var lightDirection = lightPosition - point;

            normal.Normalize();
            lightDirection.Normalize();

            return OverZero(Vector3.Dot(normal, lightDirection));
        }

        public Vector3 Compute_R(Vector3 worldPosition, Vector3 normal, Vector3 lightPosition)
        {
            var L = lightPosition - worldPosition;
            L.Normalize();
            var N = normal;
            N.Normalize();
            var LoN = Vector3.Dot(L, N);
            return N * (2 * LoN) - L;
        }

        public Vector3 Compute_H(Vector3 worldPosition, Vector3 viewerPosition, Vector3 lightPosition)
        {
            var L = lightPosition - worldPosition;
            var V = viewerPosition - worldPosition;
            var SumLV = L + V;
            var len = SumLV.Length;
            var H = SumLV / len;
            return H;
        }

        #region FlatShading
        void ProcessScanLine_FS(ScanLineData data, Vector3 A, Vector3 B, Vector3 C, Vector3 D, Color color, bool blin)
        {
            var gradient1 = A.y != B.y ? (data.currentY - A.y) / (B.y - A.y) : 1;
            var gradient2 = C.y != D.y ? (data.currentY - C.y) / (D.y - C.y) : 1;

            int sx = (int)Interpolate(A.x, B.x, gradient1);
            int ex = (int)Interpolate(C.x, D.x, gradient2);

            double z1 = Interpolate(A.z, B.z, gradient1);
            double z2 = Interpolate(C.z, D.z, gradient2);

            for (var x = sx; x < ex; x++)
            {
                double gradient = (x - sx) / (double)(ex - sx);
                var z = Interpolate(z1, z2, gradient);
                var NoL = data.NoLa;
                var RoV = data.RoVa;
                var NoH = data.NoHa;
                byte r = 0, g = 0, b = 0;

                if(blin)
                {
                    r = Limit(kd * color.R * NoL + ks * color.R * NoH);
                    g = Limit(kd * color.R * NoL + ks * color.R * NoH);
                    b = Limit(kd * color.R * NoL + ks * color.R * NoH);
                }
                else
                {
                    r = Limit(kd * color.R * NoL + ks * color.R * RoV);
                    g = Limit(kd * color.R * NoL + ks * color.R * RoV);
                    b = Limit(kd * color.R * NoL + ks * color.R * RoV);
                }

                var shadedColor = Color.FromArgb(color.A, r, g, b);

                DrawPoint(new Vector3(x, data.currentY, z), shadedColor);
            }
        }

        public void DrawTriangle_FS(Vertex v1, Vertex v2, Vertex v3, Color color, double Time, Vector3 viewerPosition, bool blin = false)
        {
            if (v1.Position.y > v2.Position.y)
                Swap(ref v1, ref v2);

            if (v2.Position.y > v3.Position.y)
                Swap(ref v2, ref v3);

            if (v1.Position.y > v2.Position.y)
                Swap(ref v1, ref v2);

            var p1 = v1.Position;
            var p2 = v2.Position;
            var p3 = v3.Position;

            Vector3 vnFace = (v1.Normal + v2.Normal + v3.Normal) / 3;
            Vector3 centerPoint = (v1.WorldPosition + v2.WorldPosition + v3.WorldPosition) / 3;
            Vector3 lightPos = LP;

            double NoL = Compute_NoL(centerPoint, vnFace, lightPos);

            double NoH = 0, RoV = 0;
            if(blin)
            {
                var H = Compute_H(v1.WorldPosition, viewerPosition, lightPos);
                H.Normalize();
                var Nnorm1 = v1.Normal;
                Nnorm1.Normalize();
                NoH = Vector3.Dot(Nnorm1, H);
                NoH = Math.Pow(NoH, alpha);
            }
            else
            {
                var R = Compute_R(centerPoint, vnFace, lightPos);
                R.Normalize();
                var V = viewerPosition - centerPoint;
                V.Normalize();
                RoV = OverZero(Vector3.Dot(R, V));
                RoV = Math.Pow(RoV, alpha);
            }

            var data = new ScanLineData { NoLa = NoL, RoVa = RoV, NoHa = NoH };

            double dP1P2, dP1P3;

            if (p2.y - p1.y > 0)
                dP1P2 = (p2.x - p1.x) / (p2.y - p1.y);
            else
                dP1P2 = 0;

            if (p3.y - p1.y > 0)
                dP1P3 = (p3.x - p1.x) / (p3.y - p1.y);
            else
                dP1P3 = 0;

            if (dP1P2 > dP1P3)
                for (var y = (int)p1.y; y <= (int)p3.y; y++)
                {
                    data.currentY = y;
                    if (y < p2.y)
                        ProcessScanLine_FS(data, p1, p3, p1, p2, color, blin);
                    else
                        ProcessScanLine_FS(data, p1, p3, p2, p3, color, blin);
                }
            else
                for (var y = (int)p1.y; y <= (int)p3.y; y++)
                {
                    data.currentY = y;
                    if (y < p2.y)
                        ProcessScanLine_FS(data, p1, p2, p1, p3, color, blin);
                    else
                        ProcessScanLine_FS(data, p2, p3, p1, p3, color, blin);
                }

        }

        #endregion

        #region GoraudShading
        void ProcessScanLine_GS(ScanLineData data, Vertex A, Vertex B, Vertex C, Vertex D, Color color, bool blin)
        {
            var pa = A.Position;
            var pb = B.Position;
            var pc = C.Position;
            var pd = D.Position;

            var gradient1 = pa.y != pb.y ? (data.currentY - pa.y) / (pb.y - pa.y) : 1;
            var gradient2 = pc.y != pd.y ? (data.currentY - pc.y) / (pd.y - pc.y) : 1;

            int sx = (int)Interpolate(pa.x, pb.x, gradient1);
            int ex = (int)Interpolate(pc.x, pd.x, gradient2);

            double z1 = Interpolate(pa.z, pb.z, gradient1);
            double z2 = Interpolate(pc.z, pd.z, gradient2);

            var sNoL = Interpolate(data.NoLa, data.NoLb, gradient1);
            var eNoL = Interpolate(data.NoLc, data.NoLd, gradient2);

            double sNoH = 0, eNoH = 0, sRoV = 0, eRoV = 0;
            if (blin)
            {
                sNoH = Interpolate(data.NoHa, data.NoHb, gradient1);
                eNoH = Interpolate(data.NoHc, data.NoHd, gradient2);
            }
            else
            {
                sRoV = Interpolate(data.RoVa, data.RoVb, gradient1);
                eRoV = Interpolate(data.RoVc, data.RoVd, gradient2);
            }


            for (var x = sx; x < ex; x++)
            {
                double gradient = (x - sx) / (double)(ex - sx);
                var z = Interpolate(z1, z2, gradient);
                var NoL = Interpolate(sNoL, eNoL, gradient);
                byte r = 0, g = 0, b = 0;
                if (blin)
                {
                    var NoH = Interpolate(sNoH, eNoH, gradient);
                    r = Limit(kd * color.R * NoL + ks * color.R * NoH);
                    g = Limit(kd * color.R * NoL + ks * color.R * NoH);
                    b = Limit(kd * color.R * NoL + ks * color.R * NoH);
                }
                else
                {
                    var RoV = Interpolate(sRoV, eRoV, gradient);
                    r = Limit(kd * color.R * NoL + ks * color.R * RoV);
                    g = Limit(kd * color.R * NoL + ks * color.R * RoV);
                    b = Limit(kd * color.R * NoL + ks * color.R * RoV);
                }

                var shadedColor = Color.FromArgb(color.A, r,g,b);
                DrawPoint(new Vector3(x, data.currentY, z), shadedColor);
            }
        }

        public void DrawTriangle_GS(Vertex v1, Vertex v2, Vertex v3, Color color, double Time, Vector3 viewerPosition, bool blin = false)
        {
            if (v1.Position.y > v2.Position.y)
                Swap(ref v1, ref v2);

            if (v2.Position.y > v3.Position.y)
                Swap(ref v2, ref v3);

            if (v1.Position.y > v2.Position.y)
                Swap(ref v1, ref v2);

            var p1 = v1.Position;
            var p2 = v2.Position;
            var p3 = v3.Position;

            Vector3 lightPos = LP;

            double NoL1 = Compute_NoL(v1.WorldPosition, v1.Normal, lightPos);
            double NoL2 = Compute_NoL(v2.WorldPosition, v2.Normal, lightPos);
            double NoL3 = Compute_NoL(v3.WorldPosition, v3.Normal, lightPos);

            double RoV1 = 0, RoV2 = 0, RoV3 = 0;
            double NoH1 = 0, NoH2 = 0, NoH3 = 0;

            if (blin)
            {
                var H1 = Compute_H(v1.WorldPosition, viewerPosition, lightPos);
                var H2 = Compute_H(v2.WorldPosition, viewerPosition, lightPos);
                var H3 = Compute_H(v3.WorldPosition, viewerPosition, lightPos);

                H1.Normalize();
                H2.Normalize();
                H3.Normalize();

                var Nnorm1 = v1.Normal;
                var Nnorm2 = v2.Normal;
                var Nnorm3 = v3.Normal;

                Nnorm1.Normalize();
                Nnorm2.Normalize();
                Nnorm3.Normalize();

                NoH1 = Vector3.Dot(Nnorm1, H1);
                NoH2 = Vector3.Dot(Nnorm2, H2);
                NoH3 = Vector3.Dot(Nnorm3, H3);

                NoH1 = Math.Pow(NoH1, alpha);
                NoH2 = Math.Pow(NoH2, alpha);
                NoH3 = Math.Pow(NoH3, alpha);
            }
            else
            {
                var R1 = Compute_R(v1.WorldPosition, v1.Normal, lightPos);
                var R2 = Compute_R(v2.WorldPosition, v2.Normal, lightPos);
                var R3 = Compute_R(v3.WorldPosition, v3.Normal, lightPos);
                R1.Normalize();
                R2.Normalize();
                R3.Normalize();

                var V1 = viewerPosition - v1.WorldPosition;
                var V2 = viewerPosition - v2.WorldPosition;
                var V3 = viewerPosition - v3.WorldPosition;
                V1.Normalize();
                V2.Normalize();
                V3.Normalize();

                RoV1 = OverZero(Vector3.Dot(R1, V1));
                RoV2 = OverZero(Vector3.Dot(R2, V2));
                RoV3 = OverZero(Vector3.Dot(R3, V3));

                RoV1 = Math.Pow(RoV1, alpha);
                RoV2 = Math.Pow(RoV2, alpha);
                RoV3 = Math.Pow(RoV3, alpha);
            }

            var data = new ScanLineData { };

            double dP1P2, dP1P3;

            if (p2.y - p1.y > 0)
                dP1P2 = (p2.x - p1.x) / (p2.y - p1.y);
            else
                dP1P2 = 0;

            if (p3.y - p1.y > 0)
                dP1P3 = (p3.x - p1.x) / (p3.y - p1.y);
            else
                dP1P3 = 0;

            if (dP1P2 > dP1P3)
                for (var y = (int)p1.y; y <= (int)p3.y; y++)
                {
                    data.currentY = y;
                    if (y < p2.y)
                    {
                        data.NoLa = NoL1;
                        data.NoLb = NoL3;
                        data.NoLc = NoL1;
                        data.NoLd = NoL2;
                        if(blin)
                        {
                            data.NoHa = NoH1;
                            data.NoHb = NoH3;
                            data.NoHc = NoH1;
                            data.NoHd = NoH2;
                        }
                        else
                        {
                            data.RoVa = RoV1;
                            data.RoVb = RoV3;
                            data.RoVc = RoV1;
                            data.RoVd = RoV2;
                        }

                        ProcessScanLine_GS(data, v1, v3, v1, v2, color,blin);
                    }
                    else
                    {
                        data.NoLa = NoL1;
                        data.NoLb = NoL3;
                        data.NoLc = NoL2;
                        data.NoLd = NoL3;
                        if(blin)
                        {
                            data.NoHa = NoH1;
                            data.NoHb = NoH3;
                            data.NoHc = NoH2;
                            data.NoHd = NoH3;
                        }
                        else
                        {
                            data.RoVa = RoV1;
                            data.RoVb = RoV3;
                            data.RoVc = RoV2;
                            data.RoVd = RoV3;
                        }
           
                        ProcessScanLine_GS(data, v1, v3, v2, v3, color,blin);
                    }
                }
            else
                for (var y = (int)p1.y; y <= (int)p3.y; y++)
                {
                    data.currentY = y;
                    if (y < p2.y)
                    {
                        data.NoLa = NoL1;
                        data.NoLb = NoL2;
                        data.NoLc = NoL1;
                        data.NoLd = NoL3;
                        if(blin)
                        {
                            data.NoHa = NoH1;
                            data.NoHb = NoH2;
                            data.NoHc = NoH1;
                            data.NoHd = NoH3;
                        }
                        else
                        {
                            data.RoVa = RoV1;
                            data.RoVb = RoV2;
                            data.RoVc = RoV1;
                            data.RoVd = RoV3;
                        }
                      
                        ProcessScanLine_GS(data, v1, v2, v1, v3, color,blin);
                    }
                    else
                    {
                        data.NoLa = NoL2;
                        data.NoLb = NoL3;
                        data.NoLc = NoL1;
                        data.NoLd = NoL3;
                        if(blin)
                        {
                            data.NoHa = NoH1;
                            data.NoHb = NoH2;
                            data.NoHc = NoH1;
                            data.NoHd = NoH3;
                        }
                        else
                        {
                            data.RoVa = RoV1;
                            data.RoVb = RoV2;
                            data.RoVc = RoV1;
                            data.RoVd = RoV3;
                        }
                  
                        ProcessScanLine_GS(data, v2, v3, v1, v3, color,blin);
                    }
                }
        }
        #endregion

        #region PhongShading

        void ProcessScanLine_PS(ScanLineDataPhong data, Vertex A, Vertex B, Vertex C, Vertex D, Color color, Vector3 lightPosition, Vector3 viewerPosition, bool blin)
        {
            var pa = A.Position;
            var pb = B.Position;
            var pc = C.Position;
            var pd = D.Position;

            var wpa = A.WorldPosition;
            var wpb = B.WorldPosition;
            var wpc = C.WorldPosition;
            var wpd = D.WorldPosition;

            var gradient1 = pa.y != pb.y ? (data.currentY - pa.y) / (pb.y - pa.y) : 1;
            var gradient2 = pc.y != pd.y ? (data.currentY - pc.y) / (pd.y - pc.y) : 1;

            int sx = (int)Interpolate(pa.x, pb.x, gradient1);
            int ex = (int)Interpolate(pc.x, pd.x, gradient2);
            int swx = (int)Interpolate(wpa.x, wpb.x, gradient1);
            int ewx = (int)Interpolate(wpc.x, wpd.x, gradient2);

            double z1 = Interpolate(pa.z, pb.z, gradient1);
            double z2 = Interpolate(pc.z, pd.z, gradient2);
            double wz1 = Interpolate(wpa.z, wpb.z, gradient1);
            double wz2 = Interpolate(wpc.z, wpd.z, gradient2);

            var sN = InterpolateVector(data.Na, data.Nb, gradient1);
            var eN = InterpolateVector(data.Nc, data.Nd, gradient2);

            for (var x = sx; x < ex; x++)
            {
                double gradient = (x - sx) / (double)(ex - sx);
                var z = Interpolate(z1, z2, gradient);
                var position = new Vector3(x, data.currentY, z);

                var wx = Interpolate(swx, ewx, gradient);
                var wy = data.currentWorldY;
                var wz = Interpolate(wz1, wz2, gradient);
                var wPosition = new Vector3(wx, wy, wz);

                var N = InterpolateVector(sN, eN, gradient);
                var NoL = Compute_NoL(wPosition, N, lightPosition);

                byte r = 0, g = 0, b = 0;
                if (blin)
                {
                    var H = Compute_H(wPosition, viewerPosition, lightPosition);
                    H.Normalize();
                    var Nnorm = N;
                    Nnorm.Normalize();
                    var NoH = Vector3.Dot(Nnorm, H);
                    NoH = Math.Pow(NoH, alpha);
                    r = Limit(kd * color.R * NoL + ks * color.R * NoH);
                    g = Limit(kd * color.R * NoL + ks * color.R * NoH);
                    b = Limit(kd * color.R * NoL + ks * color.R * NoH);
                }
                else
                {
                    var R = Compute_R(wPosition, N, lightPosition);
                    R.Normalize();
                    var V = viewerPosition - wPosition;
                    V.Normalize();
                    var RoV = OverZero(Vector3.Dot(R, V));
                    RoV = Math.Pow(RoV, alpha);
                    r = Limit(kd * color.R * NoL + ks * color.R * RoV);
                    g = Limit(kd * color.R * NoL + ks * color.R * RoV);
                    b = Limit(kd * color.R * NoL + ks * color.R * RoV);
                }

             

                var shadedColor = Color.FromArgb(color.A, r, g, b);
                DrawPoint(position, shadedColor);
            }
        }

        public void DrawTriangle_PS(Vertex v1, Vertex v2, Vertex v3, Color color, double Time, Vector3 viewerPos, bool blin = false)
        {
            if (v1.Position.y > v2.Position.y)
                Swap(ref v1, ref v2);

            if (v2.Position.y > v3.Position.y)
                Swap(ref v2, ref v3);

            if (v1.Position.y > v2.Position.y)
                Swap(ref v1, ref v2);

            var p1 = v1.Position;
            var p2 = v2.Position;
            var p3 = v3.Position;

            var wp1 = v1.WorldPosition;
            var wp2 = v2.WorldPosition;
            var wp3 = v3.WorldPosition;

            Vector3 lightPos = LP;

            var data = new ScanLineDataPhong { };

            double dP1P2, dP1P3;

            if (p2.y - p1.y > 0)
                dP1P2 = (p2.x - p1.x) / (p2.y - p1.y);
            else
                dP1P2 = 0;

            if (p3.y - p1.y > 0)
                dP1P3 = (p3.x - p1.x) / (p3.y - p1.y);
            else
                dP1P3 = 0;

            if (dP1P2 > dP1P3)
                for (var y = (int)p1.y; y <= (int)p3.y; y++)
                {
                    double gradient = (y - p1.y) / (double)(p3.y - p1.y);
                    data.currentWorldY = Interpolate(wp1.y, wp3.y, gradient);
                    data.currentY = y;
                    if (y < p2.y)
                    {
                        data.Na = v1.Normal;
                        data.Nb = v3.Normal;
                        data.Nc = v1.Normal;
                        data.Nd = v2.Normal;
                        ProcessScanLine_PS(data, v1, v3, v1, v2, color, lightPos, viewerPos,blin);
                    }
                    else
                    {
                        data.Na = v1.Normal;
                        data.Nb = v3.Normal;
                        data.Nc = v2.Normal;
                        data.Nd = v3.Normal;
                        ProcessScanLine_PS(data, v1, v3, v2, v3, color, lightPos, viewerPos,blin);
                    }
                }
            else
                for (var y = (int)p1.y; y <= (int)p3.y; y++)
                {
                    double gradient = (y - p1.y) / (double)(p3.y - p1.y);
                    data.currentWorldY = Interpolate(wp1.y, wp3.y, gradient);
                    data.currentY = y;
                    if (y < p2.y)
                    {
                        data.Na = v1.Normal;
                        data.Nb = v2.Normal;
                        data.Nc = v1.Normal;
                        data.Nd = v3.Normal;
                        ProcessScanLine_PS(data, v1, v2, v1, v3, color, lightPos, viewerPos, blin);
                    }
                    else
                    {
                        data.Na = v2.Normal;
                        data.Nb = v3.Normal;
                        data.Nc = v1.Normal;
                        data.Nd = v3.Normal;
                        ProcessScanLine_PS(data, v2, v3, v1, v3, color, lightPos, viewerPos,blin);
                    }
                }
        }
        #endregion

    }
}

