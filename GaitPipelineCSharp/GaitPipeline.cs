using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Globalization;

namespace GaitPipelineCSharp
{
    public class GaitPipeline
    {
        public const double FS = 100.0;
        public const double L_LEG = 0.45;
        public const double G = 9.81;
        public const double MIN_PROMINENCE_THETA = 0.12;
        public const int MIN_DISTANCE_SAMPLES = 60;
        public const int SG_WINDOW = 31;

        private static double TimeToSeconds(string ts)
        {
            try
            {
                string[] parts = ts.Split(':');
                int h = int.Parse(parts[0]);
                int m = int.Parse(parts[1]);
                double s = ParseEuropeanFloat(parts[2]);
                return h * 3600.0 + m * 60.0 + s;
            }
            catch
            {
                return double.NaN;
            }
        }

        public static double ParseEuropeanFloat(string s)
        {
            if (string.IsNullOrWhiteSpace(s)) return double.NaN;
            s = s.Replace(",", ".");
            if (double.TryParse(s, NumberStyles.Any, CultureInfo.InvariantCulture, out double result))
                return result;
            return double.NaN;
        }

        public static List<SensorDataRow> LoadSensorData(string filepath, string sensorCanal)
        {
            Console.WriteLine($"  Cargando {Path.GetFileName(filepath)} -> canal '{sensorCanal}'...");
            
            var rowsId1 = new List<SensorDataRow>();
            var rowsId2 = new List<SensorDataRow>();
            
            string[] lines = File.ReadAllLines(filepath);
            for (int i = 3; i < lines.Length; i++)
            {
                string line = lines[i].Trim();
                if (string.IsNullOrEmpty(line)) continue;
                string[] parts = line.Split(';');
                if (parts.Length < 19) continue;

                string canal = parts[2].Trim();
                if (canal != sensorCanal) continue;

                string idVal = parts[3].Trim();
                if (idVal != "1" && idVal != "2") continue;

                try
                {
                    var row = new SensorDataRow
                    {
                        TimeStr = parts[0].Trim(),
                        Sample = int.Parse(parts[4]),
                        Id = int.Parse(idVal),
                        Gap = ParseEuropeanFloat(parts[5]),
                        AccX = ParseEuropeanFloat(parts[6]),
                        AccY = ParseEuropeanFloat(parts[7]),
                        AccZ = ParseEuropeanFloat(parts[8]),
                        RawAccX = ParseEuropeanFloat(parts[9]),
                        RawAccY = ParseEuropeanFloat(parts[10]),
                        RawAccZ = ParseEuropeanFloat(parts[11]),
                        RawGirX = ParseEuropeanFloat(parts[12]),
                        RawGirY = ParseEuropeanFloat(parts[13]),
                        RawGirZ = ParseEuropeanFloat(parts[14]),
                        QX = ParseEuropeanFloat(parts[15]),
                        QY = ParseEuropeanFloat(parts[16]),
                        QZ = ParseEuropeanFloat(parts[17]),
                        QW = ParseEuropeanFloat(parts[18])
                    };

                    if (row.Id == 1) rowsId1.Add(row);
                    else if (row.Id == 2) rowsId2.Add(row);
                }
                catch { }
            }

            var merged = new List<SensorDataRow>();
            var dictId2 = rowsId2.ToDictionary(r => r.Sample);
            
            foreach (var qRow in rowsId1)
            {
                if (dictId2.TryGetValue(qRow.Sample, out var aRow))
                {
                    qRow.AccX = aRow.AccX;
                    qRow.AccY = aRow.AccY;
                    qRow.AccZ = aRow.AccZ;
                    qRow.RawAccX = aRow.RawAccX;
                    qRow.RawAccY = aRow.RawAccY;
                    qRow.RawAccZ = aRow.RawAccZ;
                    merged.Add(qRow);
                }
            }

            merged = merged.OrderBy(r => r.Sample).ToList();
            if (merged.Count == 0) throw new Exception($"Sin datos para sensor '{sensorCanal}'");

            double t0 = TimeToSeconds(merged[0].TimeStr);
            double tEndRaw = TimeToSeconds(merged[merged.Count - 1].TimeStr);
            
            int idx0 = 0;
            while (double.IsNaN(t0) && idx0 < merged.Count) t0 = TimeToSeconds(merged[idx0++].TimeStr);
            
            int idxE = merged.Count - 1;
            while (double.IsNaN(tEndRaw) && idxE >= 0) tEndRaw = TimeToSeconds(merged[idxE--].TimeStr);

            double tEnd = tEndRaw - t0;
            if (tEnd < -3600) tEnd += 86400; // Corrección salto de medianoche

            // Normalize quaternions and build time
            for (int i = 0; i < merged.Count; i++)
            {
                merged[i].TimeS = (merged.Count > 1) ? (tEnd * i / (merged.Count - 1)) : 0;
                double qNorm = Math.Sqrt(merged[i].QX*merged[i].QX + merged[i].QY*merged[i].QY + merged[i].QZ*merged[i].QZ + merged[i].QW*merged[i].QW);
                if (qNorm > 0)
                {
                    merged[i].QX /= qNorm;
                    merged[i].QY /= qNorm;
                    merged[i].QZ /= qNorm;
                    merged[i].QW /= qNorm;
                }
            }
            Console.WriteLine($"    -> {merged.Count} muestras cargadas (duración: {merged.Last().TimeS:F1} s)");
            return merged;
        }

        public static void FilterSignals(List<SensorDataRow> data)
        {
            int n = data.Count;
            double[] accX = new double[n]; double[] accY = new double[n]; double[] accZ = new double[n];
            double[] girX = new double[n]; double[] girY = new double[n]; double[] girZ = new double[n];
            double[] qX = new double[n]; double[] qY = new double[n]; double[] qZ = new double[n]; double[] qW = new double[n];

            for (int i = 0; i < n; i++)
            {
                accX[i] = data[i].AccX; accY[i] = data[i].AccY; accZ[i] = data[i].AccZ;
                girX[i] = data[i].RawGirX; girY[i] = data[i].RawGirY; girZ[i] = data[i].RawGirZ;
                qX[i] = data[i].QX; qY[i] = data[i].QY; qZ[i] = data[i].QZ; qW[i] = data[i].QW;
            }

            accX = MathUtils.ButterworthFiltFilt6Hz(accX);
            accY = MathUtils.ButterworthFiltFilt6Hz(accY);
            accZ = MathUtils.ButterworthFiltFilt6Hz(accZ);
            girX = MathUtils.ButterworthFiltFilt6Hz(girX);
            girY = MathUtils.ButterworthFiltFilt6Hz(girY);
            girZ = MathUtils.ButterworthFiltFilt6Hz(girZ);
            qX = MathUtils.ButterworthFiltFilt6Hz(qX);
            qY = MathUtils.ButterworthFiltFilt6Hz(qY);
            qZ = MathUtils.ButterworthFiltFilt6Hz(qZ);
            qW = MathUtils.ButterworthFiltFilt6Hz(qW);

            for (int i = 0; i < n; i++)
            {
                data[i].AccX_Filt = accX[i]; data[i].AccY_Filt = accY[i]; data[i].AccZ_Filt = accZ[i];
                data[i].RawGirX_Filt = girX[i]; data[i].RawGirY_Filt = girY[i]; data[i].RawGirZ_Filt = girZ[i];
                
                double qNorm = Math.Sqrt(qX[i]*qX[i] + qY[i]*qY[i] + qZ[i]*qZ[i] + qW[i]*qW[i]);
                if (qNorm == 0) qNorm = 1;
                data[i].QX_Filt = qX[i] / qNorm;
                data[i].QY_Filt = qY[i] / qNorm;
                data[i].QZ_Filt = qZ[i] / qNorm;
                data[i].QW_Filt = qW[i] / qNorm;
            }
        }

        public static void QuaternionToEulerXyz(List<SensorDataRow> data)
        {
            int n = data.Count;
            double[] pitchRaw = new double[n];

            for (int i = 0; i < n; i++)
            {
                MathUtils.QuaternionToEulerZYX(data[i].QX_Filt, data[i].QY_Filt, data[i].QZ_Filt, data[i].QW_Filt, out double yaw, out double pitch, out double roll);
                data[i].Yaw = yaw;
                data[i].Pitch = pitch;
                data[i].Roll = roll;
                data[i].ThetaRaw = pitch;
                pitchRaw[i] = pitch;
            }

            double[] pitchLp = MathUtils.ButterworthFiltFilt3Hz(pitchRaw);
            double[] pitchSg = MathUtils.SavitzkyGolayFilter(pitchLp, SG_WINDOW, 3);

            for (int i = 0; i < n; i++)
            {
                data[i].Theta = pitchSg[i];
                data[i].ThetaDeg = pitchSg[i] * 180.0 / Math.PI;
            }
        }

        public static GaitEventsResult DetectGaitEvents(List<SensorDataRow> data, bool isLeftLeg)
        {
            double[] theta = data.Select(r => r.Theta).ToArray();
            double[] thetaSmooth = MathUtils.SavitzkyGolayFilter(theta, 7, 3);

            var (peaksPos, _) = MathUtils.FindPeaks(thetaSmooth, MIN_PROMINENCE_THETA, MIN_DISTANCE_SAMPLES);
            
            double[] negTheta = thetaSmooth.Select(x => -x).ToArray();
            var (peaksNeg, _) = MathUtils.FindPeaks(negTheta, MIN_PROMINENCE_THETA, MIN_DISTANCE_SAMPLES);

            return new GaitEventsResult
            {
                HsIndices = peaksPos,
                ToIndices = peaksNeg,
                ThetaSmooth = thetaSmooth
            };
        }

        public static GaitEventsResult ValidateGaitEvents(GaitEventsResult events, List<SensorDataRow> data)
        {
            int[] hs = events.HsIndices;
            int[] to = events.ToIndices;
            double[] theta = events.ThetaSmooth;
            double[] t = data.Select(r => r.TimeS).ToArray();

            var validStrides = new List<GaitEvent>();

            for (int i = 0; i < hs.Length - 1; i++)
            {
                int hsStart = hs[i];
                int hsEnd = hs[i + 1];

                var toBetween = to.Where(x => x > hsStart && x < hsEnd).ToArray();
                if (toBetween.Length == 0) continue;

                int toMid = toBetween[0];

                double dtStride = t[hsEnd] - t[hsStart];
                double dtStance = t[toMid] - t[hsStart];
                double dtSwing = t[hsEnd] - t[toMid];

                if (dtStride < 0.3 || dtStride > 3.0) continue;
                if (dtStance <= 0 || dtSwing <= 0) continue;

                validStrides.Add(new GaitEvent
                {
                    HsStartIdx = hsStart,
                    ToIdx = toMid,
                    HsEndIdx = hsEnd,
                    THsStart = t[hsStart],
                    TTo = t[toMid],
                    THsEnd = t[hsEnd],
                    DtStride = dtStride,
                    DtStance = dtStance,
                    DtSwing = dtSwing,
                    ThetaHsStart = theta[hsStart],
                    ThetaTo = theta[toMid],
                    ThetaHsEnd = theta[hsEnd]
                });
            }

            int[] hsValid = validStrides.Select(s => s.HsStartIdx).Concat(validStrides.Select(s => s.HsEndIdx)).Distinct().OrderBy(x => x).ToArray();
            int[] toValid = validStrides.Select(s => s.ToIdx).ToArray();

            return new GaitEventsResult
            {
                HsIndices = hsValid,
                ToIndices = toValid,
                ThetaSmooth = theta,
                ValidStrides = validStrides
            };
        }

        public static double ComputeStrideLengthPendulum(double dtStride, double thetaRange)
        {
            double thetaMax = thetaRange / 2.0;
            double L_leg_full = 0.85;
            double stepBase = 2.0 * L_leg_full * Math.Sin(Math.Abs(thetaMax));
            double strideBase = 2.0 * stepBase;
            double tNatural = 2 * Math.PI * Math.Sqrt(L_leg_full / G);
            double k = tNatural / Math.Max(dtStride, 0.1);
            if (k < 0.5) k = 0.5;
            if (k > 3.0) k = 3.0;
            double thighScaleFactor = 2.0;
            return strideBase * k * thighScaleFactor;
        }

        public static double ComputeStepHeight(double thetaSwingMax)
        {
            return L_LEG * (1.0 - Math.Cos(Math.Abs(thetaSwingMax)));
        }

        public static double ComputeStrideVelocityFromAcc(List<SensorDataRow> data, int startIdx, int endIdx)
        {
            if (endIdx <= startIdx + 2) return double.NaN;
            
            int len = endIdx - startIdx;
            double[] dt = new double[len - 1];
            double[] ax = new double[len];

            for (int i = 0; i < len; i++)
            {
                int idx = startIdx + i;
                if (i > 0) dt[i - 1] = data[idx].TimeS - data[idx - 1].TimeS;

                var global = MathUtils.RotateVector(data[idx].QX_Filt, data[idx].QY_Filt, data[idx].QZ_Filt, data[idx].QW_Filt, data[idx].AccX_Filt, data[idx].AccY_Filt, data[idx].AccZ_Filt);
                ax[i] = global.x;
            }

            double[] vel = new double[len];
            for (int i = 1; i < len; i++)
            {
                vel[i] = vel[i - 1] + ax[i - 1] * dt[i - 1];
            }

            double driftSlope = (vel[len - 1] - vel[0]) / (len - 1);
            double sumVel = 0;
            for (int i = 0; i < len; i++)
            {
                double corrected = vel[i] - (vel[0] + driftSlope * i);
                sumVel += Math.Abs(corrected);
            }
            return sumVel / len;
        }

        public static List<StrideMetrics> ComputeStrideMetrics(List<SensorDataRow> data, List<GaitEvent> validStrides, string label)
        {
            var records = new List<StrideMetrics>();
            double[] thetaSignal = data.Select(r => r.Theta).ToArray();

            for (int i = 0; i < validStrides.Count; i++)
            {
                var stride = validStrides[i];
                double thetaRange = Math.Abs(stride.ThetaTo - stride.ThetaHsStart);

                double thetaBaseline = stride.ThetaHsStart;
                int toIdx = stride.ToIdx;
                int hsIdx = stride.HsEndIdx;
                double thetaSwingMax = 0;

                if (hsIdx > toIdx + 1)
                {
                    double maxDiff = 0;
                    for (int j = toIdx; j <= hsIdx; j++)
                    {
                        double diff = Math.Abs(thetaSignal[j] - thetaBaseline);
                        if (diff > maxDiff) maxDiff = diff;
                    }
                    thetaSwingMax = maxDiff;
                }
                else
                {
                    thetaSwingMax = Math.Abs(stride.ThetaTo - thetaBaseline);
                }

                double stepHeight = ComputeStepHeight(thetaSwingMax);
                double strideLength = ComputeStrideLengthPendulum(stride.DtStride, thetaRange);
                double strideVelKinematic = stride.DtStride > 0 ? strideLength / stride.DtStride : double.NaN;

                records.Add(new StrideMetrics
                {
                    StrideId = i + 1,
                    Sensor = label,
                    TStart = stride.THsStart,
                    TEnd = stride.THsEnd,
                    DtStrideS = stride.DtStride,
                    DtStanceS = stride.DtStance,
                    DtSwingS = stride.DtSwing,
                    StancePct = 100.0 * stride.DtStance / stride.DtStride,
                    CadenceStepMin = 60.0 / stride.DtStride,
                    ThetaRangeRad = thetaRange,
                    ThetaRangeDeg = thetaRange * 180 / Math.PI,
                    ThetaSwingMaxDeg = thetaSwingMax * 180 / Math.PI,
                    StrideLengthM = strideLength,
                    StepHeightM = stepHeight,
                    StepHeightCm = stepHeight * 100,
                    StrideVelKinMs = strideVelKinematic
                });
            }
            return records;
        }

        public static List<StrideMetrics> EnforceAlternatingSteps(List<StrideMetrics> allSteps)
        {
            if (allSteps.Count == 0) return allSteps;

            var df = allSteps.OrderBy(s => s.TStart).ToList();
            var kept = new List<StrideMetrics>();
            string lastSensor = null;

            int i = 0;
            while (i < df.Count)
            {
                string currentSensor = df[i].Sensor;
                if (lastSensor == null || currentSensor != lastSensor)
                {
                    kept.Add(df[i]);
                    lastSensor = currentSensor;
                    i++;
                }
                else
                {
                    int j = i;
                    var conflictGroup = new List<StrideMetrics> { kept.Last() };
                    while (j < df.Count && df[j].Sensor == currentSensor)
                    {
                        conflictGroup.Add(df[j]);
                        j++;
                    }

                    var bestInGroup = conflictGroup.OrderByDescending(s => s.ThetaRangeRad).First();
                    kept[kept.Count - 1] = bestInGroup;
                    lastSensor = bestInGroup.Sensor;
                    i = j;
                }
            }
            return kept;
        }

        public static List<StrideMetrics> ReconstructGlobalTrajectory(List<StrideMetrics> stridesRight, List<StrideMetrics> stridesLeft)
        {
            var allSteps = stridesRight.Concat(stridesLeft).ToList();
            allSteps = EnforceAlternatingSteps(allSteps);

            double cumulDist = 0;
            for (int i = 0; i < allSteps.Count; i++)
            {
                cumulDist += allSteps[i].StrideLengthM / 2.0;
                allSteps[i].CumulDistanceM = cumulDist;
                allSteps[i].StepOrder = i + 1;
            }
            return allSteps;
        }

        public static GlobalMetrics ComputeGlobalMetrics(List<StrideMetrics> trajectory)
        {
            if (trajectory.Count == 0) return new GlobalMetrics();

            double totalTime = trajectory.Max(s => s.TEnd) - trajectory.Min(s => s.TStart);
            double totalDistance = trajectory.Sum(s => s.StrideLengthM / 2.0);
            double avgVelocity = totalTime > 0 ? totalDistance / totalTime : double.NaN;

            return new GlobalMetrics
            {
                TotalDistanceM = totalDistance,
                TotalTimeS = totalTime,
                AvgVelocityMs = avgVelocity,
                AvgVelocityKmh = avgVelocity * 3.6,
                NStridesTotal = trajectory.Count,
                AvgCadenceMin = trajectory.Average(s => s.CadenceStepMin),
                AvgStrideLengthM = trajectory.Average(s => s.StrideLengthM)
            };
        }

        public static (double[] t, double[] h) ComputeContinuousFootHeight(List<SensorDataRow> df, GaitEventsResult valEvents)
        {
            double[] t = df.Select(r => r.TimeS).ToArray();
            double[] theta = df.Select(r => r.Theta).ToArray();
            int[] hsIdx = valEvents.HsIndices;
            double[] h = new double[t.Length];

            if (hsIdx.Length < 2) return (t, h);

            double L_eff = L_LEG * 0.5;
            int windowMin = 101;
            double[] baselineRaw = MathUtils.MinimumFilter1d(theta, windowMin);
            
            int winBase = Math.Min(301, (theta.Length - 1) | 1);
            double[] baselineSmooth = winBase > 5 ? MathUtils.SavitzkyGolayFilter(baselineRaw, winBase, 3) : baselineRaw;

            for (int i = 0; i < theta.Length; i++)
            {
                double dtheta = theta[i] - baselineSmooth[i];
                if (dtheta < 0) dtheta = 0;
                h[i] = 2.5 * L_eff * Math.Sin(dtheta);
            }
            return (t, h);
        }

        public static double[] TimeToDistanceMap(List<StrideMetrics> trajectory, double[] t)
        {
            if (trajectory.Count == 0) return new double[t.Length];

            double totalTime = trajectory.Max(s => s.TEnd) - trajectory.Min(s => s.TStart);
            double totalDist = trajectory.Sum(s => s.StrideLengthM / 2.0);
            double avgVel = totalTime > 0 ? totalDist / totalTime : 0.0;
            double tMin = trajectory.Min(s => s.TStart);

            double[] d = new double[t.Length];
            for (int i = 0; i < t.Length; i++)
            {
                d[i] = (t[i] - tMin) * avgVel;
                if (d[i] < 0) d[i] = 0;
            }
            return d;
        }
    }
}
