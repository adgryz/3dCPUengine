using MathNet.Numerics.LinearAlgebra.Double;
using System;
using Windows.UI.Xaml;
using Windows.UI.Xaml.Controls;
using Windows.UI.Xaml.Media;
using Windows.UI.Xaml.Media.Imaging;
using static ENGINE.EnumModes;

namespace ENGINE
{

    public sealed partial class MainPage : Page
    {

        private Device device;

        private Camera camera;
        Vector3 cameraPosition;
        Vector3 cameraTargetPosition;

        Mesh car;
        Mesh donut;
        Mesh[] meshes;
        public double Time;

        CameraMode cameraMode;
        ReflectionMode reflectionMode;
        ShadingMode shadingMode;

        bool IsCarShown = true;


        public MainPage()
        {
            this.InitializeComponent();
        }

        private async void Page_Loaded(object sender, RoutedEventArgs e)
        {
            WriteableBitmap bmp = new WriteableBitmap(640, 480);
            frontBuffer.Source = bmp;

            cameraMode = CameraMode.Static;
            shadingMode = ShadingMode.Flat;
            reflectionMode = ReflectionMode.Blin;

            device = new Device(bmp, shadingMode, reflectionMode);
            SetStaticCamera();

            var carLoaded = await device.LoadJSONFileAsync("ms-appx:///Assets/car6.babylon");
            var donutLoaded = await device.LoadJSONFileAsync("ms-appx:///Assets/torus.babylon");
            car = carLoaded[0];
            donut = donutLoaded[0];
            meshes = new Mesh[] { car, donut };

            Time = 0;

            CompositionTarget.Rendering += CompositionTarget_Rendering;
        }

        private void SetStaticCamera()
        {
            cameraPosition = new Vector3(6, 0.5, 0.5);
            cameraTargetPosition = new Vector3(0, 0.5, 0.5);
            camera = new Camera(cameraPosition, cameraTargetPosition);
        }

        DateTime prevTime;
        int frameNumber = 0;
        double frameCumulative = 0;
        private void CompositionTarget_Rendering(object sender, object e)
        {
            CalculateFPS();

            device.Clear(121, 85, 72, 255);
            Time += 0.01f;

            car.ModelMatrix = GetTraversingModelMatrix(Time);
            donut.ModelMatrix = GetRotatingModelMatrix(Time);

            switch (cameraMode)
            {
                case CameraMode.Static:
                    break;
                case CameraMode.Targeting:
                    LookAtTheObject(car.ModelMatrix[0, 3], car.ModelMatrix[1, 3], car.ModelMatrix[2, 3]);
                    break;
                case CameraMode.Following:
                    FollowTheObject(car.ModelMatrix[0, 3], car.ModelMatrix[1, 3], car.ModelMatrix[2, 3], Time);
                    break;
            }

            device.Render(Time, camera, meshes);
            device.DrawScene();
        }

        private void CalculateFPS()
        {
            var time = DateTime.Now;
            var tmpFps = 1000.0 / (time - prevTime).TotalMilliseconds;
            prevTime = time;

            fps.Text = ((int)tmpFps).ToString() + " fps";

            frameNumber++;
            frameCumulative += tmpFps;
            avgFps.Text = "Current FPS " + string.Format("{0:0} fps", frameCumulative / frameNumber);
            avgFps.Text = "Average FPS " + ((int)(frameCumulative / frameNumber)).ToString() + " fps";

        }

        private DenseMatrix GetRotatingModelMatrix(double fi)
        {
            return DenseMatrix.OfArray(new double[,]
            {
                {Math.Cos(fi), -Math.Sin(fi), 0, 0},
                {Math.Sin(fi), Math.Cos(fi), 0,0 },
                {0, 0, 1, 0 },
                { 0,0,0,1 }
            });
        }
        private DenseMatrix GetTraversingModelMatrix(double fi)
        {
            int r = 2;
            return DenseMatrix.OfArray(new double[,]
            {
                {1,0, 0, -r*Math.Sin(fi)},
                {0,1, 0, -r*Math.Cos(fi) +2},
                {0, 0, 1, 0 },
                { 0,0,0,1 }
            });
        }
        private void LookAtTheObject(double objX, double objY, double objZ)
        {
            var newTarget = new Vector3(objX, objY, objZ);
            camera = new Camera(camera.Position, newTarget);
        }
        private void FollowTheObject(double objX, double objY, double objZ, double fi = 0)
        {
            var newPos = new Vector3(objX - 4, objY, objZ);
            var newTarget = new Vector3(objX - 10, objY, objZ);
            camera = new Camera(newPos, newTarget);
        }

        private void Static_Click(object sender, RoutedEventArgs e) { SetStaticCamera(); cameraMode = CameraMode.Static; Time = 0; }
        private void Targeting_Click(object sender, RoutedEventArgs e) { SetStaticCamera(); cameraMode = CameraMode.Targeting; Time = 0; }
        private void Following_Click(object sender, RoutedEventArgs e) { SetStaticCamera(); cameraMode = CameraMode.Following; Time = 0; }

        private void Flat_Click(object sender, RoutedEventArgs e) => device.shadingMode = ShadingMode.Flat;
        private void Goraud_Click(object sender, RoutedEventArgs e) => device.shadingMode = ShadingMode.Goraud;
        private void PhongSh_Click(object sender, RoutedEventArgs e) => device.shadingMode = ShadingMode.Phong;

        private void PhongSp_Click(object sender, RoutedEventArgs e) => device.reflectionMode = ReflectionMode.Phong;
        private void BlinnPhong_Click(object sender, RoutedEventArgs e) => device.reflectionMode = ReflectionMode.Blin;

        private void Car_Click(object sender, RoutedEventArgs e)
        {
            IsCarShown = !IsCarShown;

            if (IsCarShown)
                meshes = new Mesh[] { car, donut };
            else
                meshes = new Mesh[] { donut };
        }
        private void ZoomIn_Click(object sender, RoutedEventArgs e)
        {
            var newPos = new Vector3(camera.Position.x-0.5, camera.Position.y, camera.Position.z);
            camera = new Camera(newPos, camera.TargetPosition);
         }
        private void ZoomOut_Click(object sender, RoutedEventArgs e)
        {
            var newPos = new Vector3(camera.Position.x+0.5, camera.Position.y, camera.Position.z);
            camera = new Camera(newPos, camera.TargetPosition);
        }
    }
}
