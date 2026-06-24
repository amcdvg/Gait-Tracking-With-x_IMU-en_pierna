using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace GaitPipelineCSharp
{
    class Program
    {
        static void Main(string[] args)
        {
            string baseDir = "/Users/amc/Documents/gait_pierna";
            string saveDir = Path.Combine(baseDir, "resultados_csharp");
            Directory.CreateDirectory(saveDir);

            var trials = new Dictionary<string, string>
            {
                { "marcha_5kmh", "data/marcha_5kmh_50metro_jesus_20260611_181630_tem2.txt" },
                { "carrera_10kmh", "data/carrera_10kmh_100metros_jesus_20260611_181807_tem2.txt" },
                { "carrera_15kmh", "data/carrera_15kmh_100metros_jesus_20260611_181919_tem2.txt" }
            };

            var allResults = new Dictionary<string, (List<StrideMetrics> right, List<StrideMetrics> left, List<StrideMetrics> traj, GlobalMetrics global)>();

            foreach (var kvp in trials)
            {
                string trialName = kvp.Key;
                string filepath = Path.Combine(baseDir, kvp.Value);

                if (!File.Exists(filepath))
                {
                    Console.WriteLine($"[WARN] Archivo no encontrado: {filepath}");
                    continue;
                }

                try
                {
                    Console.WriteLine($"\n#################################################################");
                    Console.WriteLine($"  PROCESANDO: {trialName.ToUpper()}");
                    Console.WriteLine($"#################################################################");

                    // Right Leg (Mov)
                    string labelR = $"{trialName} | pierna derecha (Mov)";
                    Console.WriteLine($"\n--- {labelR} ---");
                    var dataR = GaitPipeline.LoadSensorData(filepath, "Mov");
                    GaitPipeline.FilterSignals(dataR);
                    GaitPipeline.QuaternionToEulerXyz(dataR);
                    var eventsR = GaitPipeline.DetectGaitEvents(dataR, false);
                    var validEventsR = GaitPipeline.ValidateGaitEvents(eventsR, dataR);
                    var metricsR = GaitPipeline.ComputeStrideMetrics(dataR, validEventsR.ValidStrides, "Mov");
                    Plotter.PlotRawVsFiltered(dataR, labelR, saveDir);
                    Plotter.PlotTheta(dataR, labelR, validEventsR, saveDir);

                    // Left Leg (Mov2)
                    string labelL = $"{trialName} | pierna izquierda (Mov2)";
                    Console.WriteLine($"\n--- {labelL} ---");
                    var dataL = GaitPipeline.LoadSensorData(filepath, "Mov2");
                    GaitPipeline.FilterSignals(dataL);
                    GaitPipeline.QuaternionToEulerXyz(dataL);
                    var eventsL = GaitPipeline.DetectGaitEvents(dataL, true);
                    var validEventsL = GaitPipeline.ValidateGaitEvents(eventsL, dataL);
                    var metricsL = GaitPipeline.ComputeStrideMetrics(dataL, validEventsL.ValidStrides, "Mov2");
                    Plotter.PlotRawVsFiltered(dataL, labelL, saveDir);
                    Plotter.PlotTheta(dataL, labelL, validEventsL, saveDir);

                    // Global Trajectory
                    var trajectory = GaitPipeline.ReconstructGlobalTrajectory(metricsR, metricsL);
                    var globalMetrics = GaitPipeline.ComputeGlobalMetrics(trajectory);

                    Console.WriteLine($"\n============================================================");
                    Console.WriteLine($"MÉTRICAS GLOBALES — {trialName.ToUpper()}");
                    Console.WriteLine($"============================================================");
                    Console.WriteLine($"  TotalDistanceM: {globalMetrics.TotalDistanceM:F3}");
                    Console.WriteLine($"  TotalTimeS: {globalMetrics.TotalTimeS:F3}");
                    Console.WriteLine($"  AvgVelocityMs: {globalMetrics.AvgVelocityMs:F3}");
                    Console.WriteLine($"  AvgVelocityKmh: {globalMetrics.AvgVelocityKmh:F3}");
                    Console.WriteLine($"  NStridesTotal: {globalMetrics.NStridesTotal}");
                    Console.WriteLine($"  AvgCadenceMin: {globalMetrics.AvgCadenceMin:F3}");
                    Console.WriteLine($"  AvgStrideLengthM: {globalMetrics.AvgStrideLengthM:F3}");

                    // Global Plots
                    Plotter.PlotCumulativeDistance(trajectory, trialName, saveDir);
                    Plotter.PlotHeightVsDistance(trajectory, trialName, saveDir);
                    Plotter.PlotContinuousHeight(dataR, dataL, validEventsR, validEventsL, trajectory, trialName, saveDir);
                    Plotter.PlotThetaComparison(dataR, dataL, validEventsR, validEventsL, trajectory, trialName, saveDir);

                    allResults[trialName] = (metricsR, metricsL, trajectory, globalMetrics);
                }
                catch (Exception ex)
                {
                    Console.WriteLine($"[ERROR] Error procesando {trialName}: {ex.Message}");
                }
            }

            if (allResults.Count > 0)
            {
                ExportResults(allResults, saveDir);
            }
            
            Console.WriteLine("\n✓ Pipeline en C# completado.");
        }

        static void ExportResults(Dictionary<string, (List<StrideMetrics> right, List<StrideMetrics> left, List<StrideMetrics> traj, GlobalMetrics global)> allResults, string outputDir)
        {
            var summaryLines = new List<string> { "trial,total_distance_m,total_time_s,avg_velocity_ms,avg_velocity_kmh,n_strides_total,avg_cadence_min,avg_stride_length_m" };

            foreach (var kvp in allResults)
            {
                string trial = kvp.Key;
                var res = kvp.Value;

                ExportStrides(res.right, Path.Combine(outputDir, $"{trial}_metrics_right.csv"), false);
                ExportStrides(res.left, Path.Combine(outputDir, $"{trial}_metrics_left.csv"), false);
                ExportStrides(res.traj, Path.Combine(outputDir, $"{trial}_trajectory.csv"), true);

                var g = res.global;
                summaryLines.Add(string.Format(System.Globalization.CultureInfo.InvariantCulture, 
                    "{0},{1},{2},{3},{4},{5},{6},{7}", 
                    trial, g.TotalDistanceM, g.TotalTimeS, g.AvgVelocityMs, g.AvgVelocityKmh, g.NStridesTotal, g.AvgCadenceMin, g.AvgStrideLengthM));
            }

            File.WriteAllLines(Path.Combine(outputDir, "global_summary.csv"), summaryLines);
            Console.WriteLine($"\n  ✓ Resultados CSV exportados a: {outputDir}");
        }

        static void ExportStrides(List<StrideMetrics> strides, string path, bool isTrajectory)
        {
            if (strides == null || strides.Count == 0) return;
            var header = "stride_id,sensor,t_start,t_end,dt_stride_s,dt_stance_s,dt_swing_s,stance_pct,cadence_step_min,theta_range_rad,theta_range_deg,theta_swing_max_deg,stride_length_m,step_height_m,step_height_cm,stride_vel_kin_ms";
            if (isTrajectory) {
                header += ",step_order,cumul_distance_m";
            }
            var lines = new List<string> { header };
            foreach (var s in strides)
            {
                string line = string.Format(System.Globalization.CultureInfo.InvariantCulture,
                    "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15}",
                    s.StrideId, s.Sensor, s.TStart, s.TEnd, s.DtStrideS, s.DtStanceS, s.DtSwingS,
                    s.StancePct, s.CadenceStepMin, s.ThetaRangeRad, s.ThetaRangeDeg, s.ThetaSwingMaxDeg,
                    s.StrideLengthM, s.StepHeightM, s.StepHeightCm, s.StrideVelKinMs);
                if (isTrajectory) {
                    line += string.Format(System.Globalization.CultureInfo.InvariantCulture, ",{0},{1}", s.StepOrder, s.CumulDistanceM);
                }
                lines.Add(line);
            }
            File.WriteAllLines(path, lines);
        }
    }
}
