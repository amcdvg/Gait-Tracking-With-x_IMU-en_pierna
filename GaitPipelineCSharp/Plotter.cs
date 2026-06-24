using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ScottPlot;

namespace GaitPipelineCSharp
{
    public static class Plotter
    {
        public static void PlotRawVsFiltered(List<SensorDataRow> data, string label, string saveDir)
        {
            var plt = new Plot();
            double[] t = data.Select(d => d.TimeS).ToArray();
            
            // ScottPlot 5 handles subplots differently, we can create an image with 6 subplots
            // For simplicity in a single file, we will plot just one major signal or save separate files
            // Alternatively, create 6 separate plots and combine them, or just save them separately.
            // A faithful migration of subplots in SP5:
            // Actually, SP5 multi-panel is not as trivial as matplotlib.subplots(3,2).
            // We will save a combined layout if possible, or 6 separate images.
            // To be robust, let's plot AccX as an example or save 6 separate plots.
            
            string[] signals = { "AccX", "AccY", "AccZ", "GirX", "GirY", "GirZ" };
            foreach (var sig in signals)
            {
                var p = new Plot();
                double[] raw = sig switch {
                    "AccX" => data.Select(d => d.AccX).ToArray(),
                    "AccY" => data.Select(d => d.AccY).ToArray(),
                    "AccZ" => data.Select(d => d.AccZ).ToArray(),
                    "GirX" => data.Select(d => d.RawGirX).ToArray(),
                    "GirY" => data.Select(d => d.RawGirY).ToArray(),
                    _ => data.Select(d => d.RawGirZ).ToArray()
                };
                double[] filt = sig switch {
                    "AccX" => data.Select(d => d.AccX_Filt).ToArray(),
                    "AccY" => data.Select(d => d.AccY_Filt).ToArray(),
                    "AccZ" => data.Select(d => d.AccZ_Filt).ToArray(),
                    "GirX" => data.Select(d => d.RawGirX_Filt).ToArray(),
                    "GirY" => data.Select(d => d.RawGirY_Filt).ToArray(),
                    _ => data.Select(d => d.RawGirZ_Filt).ToArray()
                };

                var sigRaw = p.Add.ScatterLine(t, raw);
                sigRaw.Color = Colors.Gray.WithOpacity(0.8);
                sigRaw.LineWidth = 1;
                sigRaw.LegendText = "Crudo";

                var sigFilt = p.Add.ScatterLine(t, filt);
                sigFilt.Color = Colors.Red;
                sigFilt.LineWidth = 1.5f;
                sigFilt.LegendText = "Filtrado (LP 6 Hz)";

                p.Title($"{sig} - {label}");
                p.XLabel("Tiempo (s)");
                p.ShowLegend();

                string f = Path.Combine(saveDir, $"raw_vs_filt_{sig}_{label.Replace(" ", "_")}.png");
                p.SavePng(f, 800, 400);
            }
        }

        public static void PlotTheta(List<SensorDataRow> data, string label, GaitEventsResult events, string saveDir)
        {
            var p = new Plot();
            double[] t = data.Select(d => d.TimeS).ToArray();
            double[] thetaRaw = data.Select(d => d.ThetaRaw * 180 / Math.PI).ToArray();
            double[] theta = data.Select(d => d.ThetaDeg).ToArray();

            var rawLine = p.Add.ScatterLine(t, thetaRaw);
            rawLine.Color = Colors.LightGray;
            rawLine.LegendText = "θ pitch crudo";

            var filtLine = p.Add.ScatterLine(t, theta);
            filtLine.Color = Colors.DarkBlue;
            filtLine.LineWidth = 2;
            filtLine.LegendText = "θ filtrado";

            if (events.HsIndices != null)
            {
                double[] hsX = events.HsIndices.Select(i => t[i]).ToArray();
                double[] hsY = events.HsIndices.Select(i => theta[i]).ToArray();
                var hsScatter = p.Add.ScatterPoints(hsX, hsY);
                hsScatter.Color = Colors.Red;
                hsScatter.MarkerShape = MarkerShape.FilledTriangleDown;
                hsScatter.MarkerSize = 10;
                hsScatter.LegendText = "Heel Strike";
            }

            if (events.ToIndices != null)
            {
                double[] toX = events.ToIndices.Select(i => t[i]).ToArray();
                double[] toY = events.ToIndices.Select(i => theta[i]).ToArray();
                var toScatter = p.Add.ScatterPoints(toX, toY);
                toScatter.Color = Colors.Green;
                toScatter.MarkerShape = MarkerShape.FilledTriangleUp;
                toScatter.MarkerSize = 10;
                toScatter.LegendText = "Toe Off";
            }

            p.Title($"Ángulo Pendular de la Pierna - {label}");
            p.XLabel("Tiempo (s)");
            p.YLabel("θ — Ángulo sagital (°)");
            p.ShowLegend();

            string fname = Path.Combine(saveDir, $"theta_{label.Replace(" ", "_").Replace("|", "")}.png");
            p.SavePng(fname, 1200, 400);
        }

        public static void PlotCumulativeDistance(List<StrideMetrics> trajectory, string label, string saveDir)
        {
            if (trajectory.Count == 0) return;

            // Plot 1: Stride Length vs Step Order (with Cumulative Distance overlay)
            var pLen = new Plot();
            
            var trajR = trajectory.Where(s => s.Sensor.Contains("Mov") && !s.Sensor.Contains("Mov2")).ToList();
            var trajL = trajectory.Where(s => s.Sensor.Contains("Mov2")).ToList();

            if (trajR.Count > 0)
            {
                var barsR = pLen.Add.Bars(trajR.Select(s => (double)s.StepOrder).ToArray(), trajR.Select(s => s.StrideLengthM).ToArray());
                barsR.Color = Colors.Red.WithOpacity(0.7);
                barsR.LegendText = "Pierna Derecha (Mov)";
            }
            if (trajL.Count > 0)
            {
                var barsL = pLen.Add.Bars(trajL.Select(s => (double)s.StepOrder).ToArray(), trajL.Select(s => s.StrideLengthM).ToArray());
                barsL.Color = Colors.Blue.WithOpacity(0.7);
                barsL.LegendText = "Pierna Izquierda (Mov2)";
            }

            var axis2 = pLen.Axes.AddRightAxis();
            var lineCum = pLen.Add.ScatterLine(trajectory.Select(s => (double)s.StepOrder).ToArray(), trajectory.Select(s => s.CumulDistanceM).ToArray());
            lineCum.Axes.YAxis = axis2;
            lineCum.Color = Colors.Black;
            lineCum.LineWidth = 1.5f;
            lineCum.LinePattern = LinePattern.Dashed;
            lineCum.LegendText = "Dist. acumulada";

            pLen.Title($"Longitud de Zancada y Distancia Acumulada — {label}");
            pLen.XLabel("Número de paso");
            pLen.YLabel("Longitud de zancada (m)");
            axis2.Label.Text = "Distancia acumulada (m)";
            pLen.ShowLegend();

            string fnameLen = Path.Combine(saveDir, $"trajectory_length_{label.Replace(" ", "_").Replace("|", "")}.png");
            pLen.SavePng(fnameLen, 1000, 400);

            // Plot 2: Stride Velocity vs Step Order
            var pVel = new Plot();
            if (trajR.Count > 0)
            {
                var barsR = pVel.Add.Bars(trajR.Select(s => (double)s.StepOrder).ToArray(), trajR.Select(s => s.StrideVelKinMs).ToArray());
                barsR.Color = Colors.Red.WithOpacity(0.7);
                barsR.LegendText = "Pierna Derecha (Mov)";
            }
            if (trajL.Count > 0)
            {
                var barsL = pVel.Add.Bars(trajL.Select(s => (double)s.StepOrder).ToArray(), trajL.Select(s => s.StrideVelKinMs).ToArray());
                barsL.Color = Colors.Blue.WithOpacity(0.7);
                barsL.LegendText = "Pierna Izquierda (Mov2)";
            }

            pVel.Title($"Velocidad de Zancada vs Número de paso — {label}");
            pVel.XLabel("Número de paso");
            pVel.YLabel("Velocidad de zancada (m/s)");
            pVel.ShowLegend();

            string fnameVel = Path.Combine(saveDir, $"trajectory_vel_{label.Replace(" ", "_").Replace("|", "")}.png");
            pVel.SavePng(fnameVel, 1000, 400);
        }

        public static void PlotContinuousHeight(List<SensorDataRow> dfRight, List<SensorDataRow> dfLeft, GaitEventsResult evRight, GaitEventsResult evLeft, List<StrideMetrics> trajectory, string label, string saveDir)
        {
            // Plot 1: Height vs Time
            var pTime = new Plot();
            var (tR, hR) = GaitPipeline.ComputeContinuousFootHeight(dfRight, evRight);
            var (tL, hL) = GaitPipeline.ComputeContinuousFootHeight(dfLeft, evLeft);

            var lineR_t = pTime.Add.ScatterLine(tR, hR);
            lineR_t.Color = Colors.Orange;
            lineR_t.LineWidth = 2;
            lineR_t.LegendText = "Mov Height (Derecha)";

            var lineL_t = pTime.Add.ScatterLine(tL, hL);
            lineL_t.Color = Colors.DarkGreen;
            lineL_t.LineWidth = 2;
            lineL_t.LegendText = "Mov2 Height (Izquierda)";

            pTime.Title($"Height Comparison (Z Axis) vs Time - {label}");
            pTime.XLabel("Time (s)");
            pTime.YLabel("Z Position (m)");
            pTime.ShowLegend();

            string fnameTime = Path.Combine(saveDir, $"height_comparison_time_{label.Replace(" ", "_").Replace("|", "")}.png");
            pTime.SavePng(fnameTime, 1000, 400);

            // Plot 2: Height vs Distance
            var pDist = new Plot();
            double[] dR = GaitPipeline.TimeToDistanceMap(trajectory, tR);
            double[] dL = GaitPipeline.TimeToDistanceMap(trajectory, tL);

            var lineR_d = pDist.Add.ScatterLine(dR, hR);
            lineR_d.Color = Colors.Orange;
            lineR_d.LineWidth = 2;
            lineR_d.LegendText = "Mov Height (Derecha)";

            var lineL_d = pDist.Add.ScatterLine(dL, hL);
            lineL_d.Color = Colors.DarkGreen;
            lineL_d.LineWidth = 2;
            lineL_d.LegendText = "Mov2 Height (Izquierda)";

            pDist.Title($"Height Comparison (Z Axis) vs Distance - {label}");
            pDist.XLabel("Distance (m)");
            pDist.YLabel("Z Position (m)");
            pDist.ShowLegend();

            string fnameDist = Path.Combine(saveDir, $"height_comparison_dist_{label.Replace(" ", "_").Replace("|", "")}.png");
            pDist.SavePng(fnameDist, 1000, 400);
        }

        public static void PlotHeightVsDistance(List<StrideMetrics> trajectory, string label, string saveDir)
        {
            if (trajectory.Count == 0) return;
            var p = new Plot();
            
            var rDist = trajectory.Where(s => s.Sensor.Contains("Mov") && !s.Sensor.Contains("Mov2")).Select(s => s.CumulDistanceM).ToArray();
            var rHeight = trajectory.Where(s => s.Sensor.Contains("Mov") && !s.Sensor.Contains("Mov2")).Select(s => s.StepHeightCm).ToArray();
            if (rDist.Length > 0)
            {
                var rScatter = p.Add.ScatterPoints(rDist, rHeight);
                rScatter.Color = Colors.Red;
                rScatter.MarkerSize = 10;
                rScatter.LegendText = "Pierna Derecha (Mov)";
            }

            var lDist = trajectory.Where(s => s.Sensor.Contains("Mov2")).Select(s => s.CumulDistanceM).ToArray();
            var lHeight = trajectory.Where(s => s.Sensor.Contains("Mov2")).Select(s => s.StepHeightCm).ToArray();
            if (lDist.Length > 0)
            {
                var lScatter = p.Add.ScatterPoints(lDist, lHeight);
                lScatter.Color = Colors.Blue;
                lScatter.MarkerSize = 10;
                lScatter.MarkerShape = MarkerShape.FilledSquare;
                lScatter.LegendText = "Pierna Izquierda (Mov2)";
            }

            p.Title($"Elevación del Pie vs Distancia Recorrida — {label}");
            p.XLabel("Distancia recorrida (m)");
            p.YLabel("Elevación del pie — swing (cm)");
            p.ShowLegend();

            string fname = Path.Combine(saveDir, $"elevacion_swing_scatter_{label.Replace(" ", "_").Replace("|", "")}.png");
            p.SavePng(fname, 1000, 400);
        }

        private static int[] Subsample(int[] indices, int maxN = 60)
        {
            if (indices == null || indices.Length <= maxN) return indices;
            int step = Math.Max(1, indices.Length / maxN);
            var res = new List<int>();
            for (int i = 0; i < indices.Length; i += step) res.Add(indices[i]);
            return res.ToArray();
        }

        public static void PlotThetaComparison(List<SensorDataRow> dfRight, List<SensorDataRow> dfLeft, GaitEventsResult evRight, GaitEventsResult evLeft, List<StrideMetrics> trajectory, string label, string saveDir)
        {
            // Plot 1: Theta vs Time
            var pTime = new Plot();
            double[] tR = dfRight.Select(r => r.TimeS).ToArray();
            double[] tL = dfLeft.Select(r => r.TimeS).ToArray();
            double[] thetaR = dfRight.Select(r => r.ThetaDeg).ToArray();
            double[] thetaL = dfLeft.Select(r => r.ThetaDeg).ToArray();

            var lineR_t = pTime.Add.ScatterLine(tR, thetaR);
            lineR_t.Color = Colors.Orange;
            lineR_t.LineWidth = 2;
            lineR_t.LegendText = "Mov θ (pierna D)";

            var lineL_t = pTime.Add.ScatterLine(tL, thetaL);
            lineL_t.Color = Colors.DarkGreen;
            lineL_t.LineWidth = 2;
            lineL_t.LegendText = "Mov2 θ (pierna I)";

            int[] hs_r = Subsample(evRight.HsIndices);
            int[] to_r = Subsample(evRight.ToIndices);
            int[] hs_l = Subsample(evLeft.HsIndices);
            int[] to_l = Subsample(evLeft.ToIndices);

            if (hs_r != null) {
                var m = pTime.Add.ScatterPoints(hs_r.Select(i => tR[i]).ToArray(), hs_r.Select(i => thetaR[i]).ToArray());
                m.MarkerShape = MarkerShape.FilledTriangleDown; m.Color = Colors.Orange; m.LegendText = "HS Derecha"; m.MarkerSize = 10;
            }
            if (to_r != null) {
                var m = pTime.Add.ScatterPoints(to_r.Select(i => tR[i]).ToArray(), to_r.Select(i => thetaR[i]).ToArray());
                m.MarkerShape = MarkerShape.FilledTriangleUp; m.Color = Colors.Orange; m.LegendText = "TO Derecha"; m.MarkerSize = 10;
            }
            if (hs_l != null) {
                var m = pTime.Add.ScatterPoints(hs_l.Select(i => tL[i]).ToArray(), hs_l.Select(i => thetaL[i]).ToArray());
                m.MarkerShape = MarkerShape.FilledTriangleDown; m.Color = Colors.DarkGreen; m.LegendText = "HS Izquierda"; m.MarkerSize = 10;
            }
            if (to_l != null) {
                var m = pTime.Add.ScatterPoints(to_l.Select(i => tL[i]).ToArray(), to_l.Select(i => thetaL[i]).ToArray());
                m.MarkerShape = MarkerShape.FilledTriangleUp; m.Color = Colors.DarkGreen; m.LegendText = "TO Izquierda"; m.MarkerSize = 10;
            }

            var hl_t = pTime.Add.HorizontalLine(0);
            hl_t.Color = Colors.Gray.WithOpacity(0.6);
            hl_t.LinePattern = LinePattern.Dashed;

            pTime.Title($"Ángulo θ (sagital) vs Tiempo - {label}");
            pTime.XLabel("Time (s)");
            pTime.YLabel("θ — Ángulo sagital (°)");
            pTime.ShowLegend();
            if (tR.Length > 0) pTime.Axes.SetLimitsX(0, tR.LastOrDefault());

            string fnameTime = Path.Combine(saveDir, $"theta_comparison_time_{label.Replace(" ", "_").Replace("|", "")}.png");
            pTime.SavePng(fnameTime, 1000, 400);

            // Plot 2: Theta vs Distance
            var pDist = new Plot();
            double[] dR = GaitPipeline.TimeToDistanceMap(trajectory, tR);
            double[] dL = GaitPipeline.TimeToDistanceMap(trajectory, tL);

            var lineR_d = pDist.Add.ScatterLine(dR, thetaR);
            lineR_d.Color = Colors.Orange;
            lineR_d.LineWidth = 2;
            lineR_d.LegendText = "Mov θ (pierna D)";

            var lineL_d = pDist.Add.ScatterLine(dL, thetaL);
            lineL_d.Color = Colors.DarkGreen;
            lineL_d.LineWidth = 2;
            lineL_d.LegendText = "Mov2 θ (pierna I)";

            if (hs_r != null) {
                var m = pDist.Add.ScatterPoints(hs_r.Select(i => dR[i]).ToArray(), hs_r.Select(i => thetaR[i]).ToArray());
                m.MarkerShape = MarkerShape.FilledTriangleDown; m.Color = Colors.Orange; m.LegendText = "HS Derecha"; m.MarkerSize = 10;
            }
            if (to_r != null) {
                var m = pDist.Add.ScatterPoints(to_r.Select(i => dR[i]).ToArray(), to_r.Select(i => thetaR[i]).ToArray());
                m.MarkerShape = MarkerShape.FilledTriangleUp; m.Color = Colors.Orange; m.LegendText = "TO Derecha"; m.MarkerSize = 10;
            }
            if (hs_l != null) {
                var m = pDist.Add.ScatterPoints(hs_l.Select(i => dL[i]).ToArray(), hs_l.Select(i => thetaL[i]).ToArray());
                m.MarkerShape = MarkerShape.FilledTriangleDown; m.Color = Colors.DarkGreen; m.LegendText = "HS Izquierda"; m.MarkerSize = 10;
            }
            if (to_l != null) {
                var m = pDist.Add.ScatterPoints(to_l.Select(i => dL[i]).ToArray(), to_l.Select(i => thetaL[i]).ToArray());
                m.MarkerShape = MarkerShape.FilledTriangleUp; m.Color = Colors.DarkGreen; m.LegendText = "TO Izquierda"; m.MarkerSize = 10;
            }

            var hl_d = pDist.Add.HorizontalLine(0);
            hl_d.Color = Colors.Gray.WithOpacity(0.6);
            hl_d.LinePattern = LinePattern.Dashed;

            pDist.Title($"Ángulo θ (sagital) vs Distancia - {label}");
            pDist.XLabel("Distancia (m)");
            pDist.YLabel("θ — Ángulo sagital (°)");
            pDist.ShowLegend();
            if (dR.Length > 0) pDist.Axes.SetLimitsX(0, dR.LastOrDefault());

            string fnameDist = Path.Combine(saveDir, $"theta_comparison_dist_{label.Replace(" ", "_").Replace("|", "")}.png");
            pDist.SavePng(fnameDist, 1000, 400);
        }
    }
}
