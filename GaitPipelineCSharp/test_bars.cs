using System.Collections.Generic;
using ScottPlot;
using ScottPlot.Plottables;

class Test
{
    void Run()
    {
        var p = new Plot();
        double[] x = { 1, 2, 3 };
        double[] y = { 10, 20, 30 };
        var b = p.Add.Bars(x, y);
        b.LegendText = "Test";
        
        var axis2 = p.Axes.AddRightAxis();
        var line = p.Add.Scatter(x, y);
        line.Axes.YAxis = axis2;
    }
}
