using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;

namespace GaitPipelineCSharp
{
    public static class MathUtils
    {
        // ---------------------------------------------------------------------------
        // Butterworth Filter & FiltFilt
        // ---------------------------------------------------------------------------
        
        // Coefficients derived directly from scipy.signal.butter(4, fc/(100/2))
        private static readonly double[] b_6Hz = { 0.0008063598650371013, 0.003225439460148405, 0.0048381591902226075, 0.003225439460148405, 0.0008063598650371013 };
        private static readonly double[] a_6Hz = { 1.0, -3.017555238686489, 3.5071937247162053, -1.8475509441185773, 0.3708142159294547 };
        private static readonly double[] zi_6Hz = { 0.999193640135001, -2.0215870380117513, 1.4807685275143647, -0.3700078560644317 };

        private static readonly double[] b_3Hz = { 6.23869835484794e-05, 0.0002495479341939176, 0.0003743219012908764, 0.0002495479341939176, 6.23869835484794e-05 };
        private static readonly double[] a_3Hz = { 1.0, -3.5077862073907826, 4.640902412686707, -2.7426528211203722, 0.6105348075612237 };
        private static readonly double[] zi_3Hz = { 0.9999376130169201, -2.5080981423097004, 2.13242994847789, -0.6104724205779614 };

        public static double[] ButterworthFiltFilt6Hz(double[] data)
        {
            return FiltFilt(b_6Hz, a_6Hz, zi_6Hz, data);
        }

        public static double[] ButterworthFiltFilt3Hz(double[] data)
        {
            return FiltFilt(b_3Hz, a_3Hz, zi_3Hz, data);
        }

        private static double[] LFilter(double[] b, double[] a, double[] data, double[] zi)
        {
            int n = data.Length;
            double[] y = new double[n];
            double[] z = (double[])zi.Clone();

            for (int i = 0; i < n; i++)
            {
                y[i] = b[0] * data[i] + z[0];
                for (int j = 1; j < b.Length - 1; j++)
                {
                    z[j - 1] = b[j] * data[i] + z[j] - a[j] * y[i];
                }
                z[b.Length - 2] = b[b.Length - 1] * data[i] - a[a.Length - 1] * y[i];
            }
            return y;
        }

        private static double[] FiltFilt(double[] b, double[] a, double[] zi_step, double[] data)
        {
            int n = data.Length;
            int padlen = 3 * Math.Max(a.Length, b.Length);
            
            if (n < padlen)
                throw new ArgumentException("Data too short for filtfilt");

            // Odd extension padding
            double[] ext = new double[n + 2 * padlen];
            double edgeLeft = data[0];
            double edgeRight = data[n - 1];

            for (int i = 0; i < padlen; i++)
            {
                ext[padlen - 1 - i] = 2.0 * edgeLeft - data[i + 1];
                ext[padlen + n + i] = 2.0 * edgeRight - data[n - 2 - i];
            }
            for (int i = 0; i < n; i++)
                ext[padlen + i] = data[i];

            // Initial conditions for forward filter
            double[] zi_fwd = new double[zi_step.Length];
            for (int i = 0; i < zi_step.Length; i++)
                zi_fwd[i] = zi_step[i] * ext[0];

            // Forward filter
            double[] y_fwd = LFilter(b, a, ext, zi_fwd);

            // Reverse the sequence
            Array.Reverse(y_fwd);

            // Initial conditions for backward filter
            double[] zi_bwd = new double[zi_step.Length];
            for (int i = 0; i < zi_step.Length; i++)
                zi_bwd[i] = zi_step[i] * y_fwd[0];

            // Backward filter
            double[] y_bwd = LFilter(b, a, y_fwd, zi_bwd);

            // Reverse again
            Array.Reverse(y_bwd);

            // Extract the middle part
            double[] result = new double[n];
            Array.Copy(y_bwd, padlen, result, 0, n);
            return result;
        }

        // ---------------------------------------------------------------------------
        // Savitzky-Golay Filter
        // ---------------------------------------------------------------------------
        public static double[] SavitzkyGolayFilter(double[] data, int window_length, int polyorder)
        {
            if (window_length % 2 == 0) throw new ArgumentException("window_length must be odd");
            
            int n = data.Length;
            double[] result = new double[n];
            int half_window = window_length / 2;

            // Generate SG weights for interior points
            var M = Matrix<double>.Build.Dense(window_length, polyorder + 1);
            for (int i = 0; i < window_length; i++)
            {
                double x = i - half_window;
                for (int j = 0; j <= polyorder; j++)
                    M[i, j] = Math.Pow(x, j);
            }
            var M_pinv = M.PseudoInverse();
            var weights = M_pinv.Row(0).ToArray();

            // Interior
            for (int i = half_window; i < n - half_window; i++)
            {
                double sum = 0;
                for (int j = 0; j < window_length; j++)
                {
                    sum += weights[j] * data[i - half_window + j];
                }
                result[i] = sum;
            }

            // Edges: fit polynomial to first/last window_length points
            // Left edge
            var leftX = Vector<double>.Build.Dense(window_length, i => i);
            var leftY = Vector<double>.Build.Dense(window_length, i => data[i]);
            var leftPoly = MathNet.Numerics.Fit.Polynomial(leftX.ToArray(), leftY.ToArray(), polyorder);
            for (int i = 0; i < half_window; i++)
            {
                result[i] = MathNet.Numerics.Polynomial.Evaluate(i, leftPoly);
            }

            // Right edge
            var rightX = Vector<double>.Build.Dense(window_length, i => i);
            var rightY = Vector<double>.Build.Dense(window_length, i => data[n - window_length + i]);
            var rightPoly = MathNet.Numerics.Fit.Polynomial(rightX.ToArray(), rightY.ToArray(), polyorder);
            for (int i = 0; i < half_window; i++)
            {
                int targetIdx = n - half_window + i;
                double evalX = window_length - half_window + i;
                result[targetIdx] = MathNet.Numerics.Polynomial.Evaluate(evalX, rightPoly);
            }

            return result;
        }

        // ---------------------------------------------------------------------------
        // Moving Minimum 1D (for scipy.ndimage.minimum_filter1d)
        // ---------------------------------------------------------------------------
        public static double[] MinimumFilter1d(double[] data, int size)
        {
            int n = data.Length;
            double[] result = new double[n];
            int half = size / 2;

            for (int i = 0; i < n; i++)
            {
                int start = Math.Max(0, i - half);
                int end = Math.Min(n - 1, i + half);
                double min = data[start];
                for (int j = start + 1; j <= end; j++)
                {
                    if (data[j] < min) min = data[j];
                }
                result[i] = min;
            }
            return result;
        }

        // ---------------------------------------------------------------------------
        // Peak Detection (scipy.signal.find_peaks with prominence & distance)
        // ---------------------------------------------------------------------------
        public static (int[] peaks, double[] prominences) FindPeaks(double[] data, double min_prominence, int min_distance)
        {
            int n = data.Length;
            List<int> candidate_peaks = new List<int>();

            // 1. Find local maxima
            for (int i = 1; i < n - 1; i++)
            {
                if (data[i] > data[i - 1] && data[i] > data[i + 1])
                {
                    candidate_peaks.Add(i);
                }
                // Handle flat peaks
                else if (data[i] > data[i - 1] && data[i] == data[i + 1])
                {
                    int j = i + 1;
                    while (j < n - 1 && data[j] == data[i]) j++;
                    if (data[i] > data[j])
                    {
                        candidate_peaks.Add((i + j - 1) / 2); // middle of plateau
                    }
                }
            }

            // 2. Calculate prominences
            var valid_peaks = new List<(int idx, double prom)>();
            foreach (int p in candidate_peaks)
            {
                double left_min = data[p];
                for (int i = p - 1; i >= 0; i--)
                {
                    if (data[i] > data[p]) break;
                    if (data[i] < left_min) left_min = data[i];
                }

                double right_min = data[p];
                for (int i = p + 1; i < n; i++)
                {
                    if (data[i] > data[p]) break;
                    if (data[i] < right_min) right_min = data[i];
                }

                double prom = data[p] - Math.Max(left_min, right_min);
                if (prom >= min_prominence)
                {
                    valid_peaks.Add((p, prom));
                }
            }

            // 3. Filter by distance
            if (min_distance > 1)
            {
                // Sort by prominence descending
                valid_peaks = valid_peaks.OrderByDescending(p => p.prom).ToList();
                bool[] keep = new bool[valid_peaks.Count];
                for (int i = 0; i < keep.Length; i++) keep[i] = true;

                for (int i = 0; i < valid_peaks.Count; i++)
                {
                    if (!keep[i]) continue;
                    for (int j = i + 1; j < valid_peaks.Count; j++)
                    {
                        if (keep[j] && Math.Abs(valid_peaks[i].idx - valid_peaks[j].idx) < min_distance)
                        {
                            keep[j] = false;
                        }
                    }
                }

                var filtered_peaks = new List<(int idx, double prom)>();
                for (int i = 0; i < valid_peaks.Count; i++)
                {
                    if (keep[i]) filtered_peaks.Add(valid_peaks[i]);
                }
                valid_peaks = filtered_peaks;
            }

            // Sort back by index
            valid_peaks = valid_peaks.OrderBy(p => p.idx).ToList();

            int[] final_peaks = valid_peaks.Select(p => p.idx).ToArray();
            double[] final_proms = valid_peaks.Select(p => p.prom).ToArray();
            return (final_peaks, final_proms);
        }

        // ---------------------------------------------------------------------------
        // Quaternion to Euler (ZYX) and Rotation
        // ---------------------------------------------------------------------------
        public static void QuaternionToEulerZYX(double qx, double qy, double qz, double qw, out double yaw, out double pitch, out double roll)
        {
            // scipy.spatial.transform.Rotation.as_euler("ZYX")
            // ZYX corresponds to Yaw, Pitch, Roll
            
            // Roll (x-axis rotation)
            double sinr_cosp = 2 * (qw * qx + qy * qz);
            double cosr_cosp = 1 - 2 * (qx * qx + qy * qy);
            roll = Math.Atan2(sinr_cosp, cosr_cosp);

            // Pitch (y-axis rotation)
            double sinp = 2 * (qw * qy - qz * qx);
            if (Math.Abs(sinp) >= 1)
                pitch = Math.CopySign(Math.PI / 2, sinp); // use 90 degrees if out of range
            else
                pitch = Math.Asin(sinp);

            // Yaw (z-axis rotation)
            double siny_cosp = 2 * (qw * qz + qx * qy);
            double cosy_cosp = 1 - 2 * (qy * qy + qz * qz);
            yaw = Math.Atan2(siny_cosp, cosy_cosp);
        }

        public static (double x, double y, double z) RotateVector(double qx, double qy, double qz, double qw, double vx, double vy, double vz)
        {
            // Quat multiplication: v_rot = q * v * q_conj
            // Vector part of quat q: q_v = (qx, qy, qz), scalar = qw
            // v_rot = v + 2*qw*(q_v x v) + 2*(q_v x (q_v x v))
            double tx = 2 * (qy * vz - qz * vy);
            double ty = 2 * (qz * vx - qx * vz);
            double tz = 2 * (qx * vy - qy * vx);

            double rx = vx + qw * tx + (qy * tz - qz * ty);
            double ry = vy + qw * ty + (qz * tx - qx * tz);
            double rz = vz + qw * tz + (qx * ty - qy * tx);

            return (rx, ry, rz);
        }
    }
}
