using System;
using System.Collections.Generic;

namespace GaitPipelineCSharp
{
    public class SensorDataRow
    {
        public int Sample { get; set; }
        public string TimeStr { get; set; }
        public double TimeS { get; set; }
        public int Id { get; set; }
        public double Gap { get; set; }

        public double AccX { get; set; }
        public double AccY { get; set; }
        public double AccZ { get; set; }

        public double RawAccX { get; set; }
        public double RawAccY { get; set; }
        public double RawAccZ { get; set; }

        public double RawGirX { get; set; }
        public double RawGirY { get; set; }
        public double RawGirZ { get; set; }

        public double QX { get; set; }
        public double QY { get; set; }
        public double QZ { get; set; }
        public double QW { get; set; }

        // Filtered versions
        public double AccX_Filt { get; set; }
        public double AccY_Filt { get; set; }
        public double AccZ_Filt { get; set; }
        public double RawGirX_Filt { get; set; }
        public double RawGirY_Filt { get; set; }
        public double RawGirZ_Filt { get; set; }
        public double QX_Filt { get; set; }
        public double QY_Filt { get; set; }
        public double QZ_Filt { get; set; }
        public double QW_Filt { get; set; }

        // Euler Angles and Theta
        public double Roll { get; set; }
        public double Pitch { get; set; }
        public double Yaw { get; set; }
        public double ThetaRaw { get; set; }
        public double Theta { get; set; }
        public double ThetaDeg { get; set; }
        public double ThetaSmooth { get; set; }
    }

    public class GaitEvent
    {
        public int HsStartIdx { get; set; }
        public int ToIdx { get; set; }
        public int HsEndIdx { get; set; }
        public double THsStart { get; set; }
        public double TTo { get; set; }
        public double THsEnd { get; set; }
        public double DtStride { get; set; }
        public double DtStance { get; set; }
        public double DtSwing { get; set; }
        public double ThetaHsStart { get; set; }
        public double ThetaTo { get; set; }
        public double ThetaHsEnd { get; set; }
    }

    public class GaitEventsResult
    {
        public int[] HsIndices { get; set; }
        public int[] ToIndices { get; set; }
        public double[] ThetaSmooth { get; set; }
        public List<GaitEvent> ValidStrides { get; set; }
    }

    public class StrideMetrics
    {
        public int StrideId { get; set; }
        public string Sensor { get; set; }
        public double TStart { get; set; }
        public double TEnd { get; set; }
        public double DtStrideS { get; set; }
        public double DtStanceS { get; set; }
        public double DtSwingS { get; set; }
        public double StancePct { get; set; }
        public double CadenceStepMin { get; set; }
        public double ThetaRangeRad { get; set; }
        public double ThetaRangeDeg { get; set; }
        public double ThetaSwingMaxDeg { get; set; }
        public double StrideLengthM { get; set; }
        public double StepHeightM { get; set; }
        public double StepHeightCm { get; set; }
        public double StrideVelKinMs { get; set; }
        
        // For global trajectory
        public int StepOrder { get; set; }
        public double CumulDistanceM { get; set; }
    }

    public class GlobalMetrics
    {
        public double TotalDistanceM { get; set; }
        public double TotalTimeS { get; set; }
        public double AvgVelocityMs { get; set; }
        public double AvgVelocityKmh { get; set; }
        public int NStridesTotal { get; set; }
        public double AvgCadenceMin { get; set; }
        public double AvgStrideLengthM { get; set; }
    }
}
