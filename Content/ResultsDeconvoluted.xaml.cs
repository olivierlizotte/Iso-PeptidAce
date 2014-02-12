using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
//using System.Windows.Navigation;
using System.Windows.Shapes;
using System.ComponentModel;
using System.Collections.ObjectModel;
using PeptidAce.Iso.Structures;
using FirstFloor.ModernUI.Windows;
using FirstFloor.ModernUI.Windows.Navigation;
using System.IO;

namespace PeptidAce.ModernUI.Content
{
    /// <summary>
    /// Interaction logic for About.xaml
    /// </summary>
    public partial class ResultsDeconvoluted : UserControl, IContent
    {
        private string outputFile = null;
        public static Dictionary<string, string> GetParams(string uri)
        {
            string[] splits = uri.Split(new char[] { '?', '&' });
            Dictionary<string, string> parameters = new Dictionary<string, string>();
            foreach (string str in splits)
            {
                string[] pairs = str.Split('=');
                if(pairs.Length > 1)
                    parameters.Add(pairs[0], pairs[1]);
            }
            return parameters;
        }
        public void OnFragmentNavigation(FragmentNavigationEventArgs e)
        {
        }
        public void OnNavigatedFrom(NavigationEventArgs e)
        {
        }
        public void OnNavigatedTo(NavigationEventArgs e)
        {
            StackPanel1.Children.Clear();
            Dictionary<string, string> parameters = GetParams(e.Source.OriginalString);
            if (parameters.ContainsKey("Sample") && parameters.ContainsKey("mzKey"))
            {
                foreach (Sample sample in PepIso.solver.mixedPrecursors.Keys)
                    if (sample.Name.CompareTo(parameters["Sample"]) == 0)
                        UpdateContent(sample, parameters["mzKey"]);
            }
        }
        public void OnNavigatingFrom(NavigatingCancelEventArgs e)
        {
        }

        public ResultsDeconvoluted()
        {
            InitializeComponent();
        }

        public void UpdateContent(Sample sample, string strKey)
        {    
            foreach(MixedPrecursor mixed in PepIso.solver.mixedPrecursors[sample])
            {                
                if (mixed.PeptideRatios.Count > 0 && mixed.MZ.ToString().CompareTo(strKey) == 0)
                {
                    List<double[]> intensities = new List<double[]>();
                    List<string[]> scanTimes = new List<string[]>();
                    List<double> ratios = new List<double>();
                    List<string> names = new List<string>();                    
                    //double[] timeArray = null;

                    /*
                    string[] timeStr = new string[mixed.eCurve.time.Count];

                    for (int i = 0; i < timeStr.Length; i++)
                        timeStr[i] = ((float)(mixed.eCurve.time[i] / (1000.0 * 60.0))).ToString();//ms to min
                    //scanTimes.Add(time);//

                    foreach (CharacterizedPrecursor peptide in mixed.PeptideRatios.Keys)
                    {
                        scanTimes.Add(timeStr);
                        double[] arrayInt = new double[mixed.eCurve.time.Count];
                        for (int i = 0; i < arrayInt.Length; i++)
                            //arrayInt[i] = mixed.PeptideRatios[peptide].eCurvePerMs.InterpolateIntensity(mixed.eCurve.interpolatedTime[i]);
                            arrayInt[i] = mixed.PeptideRatios[peptide].eCurveCount.InterpolateIntensity(mixed.eCurve.interpolatedTime[i]);
                        intensities.Add(arrayInt);

                        ratios.Add(mixed.PeptideRatios[peptide].eCurvePerMs.Area);

                        names.Add(peptide.Peptide.Sequence + "  [" + Utilities.vsCSV.GetFileName_NoExtension(peptide.Sample.sSDF) + "]");
                    }//*/

                    List<double> timePoints = Utilities.ElutionCurve.GetTimePoints(64, true, mixed.eCurveIntensityCount.interpolatedTime, mixed.eCurveIntensityCount.interpolatedIntensityCount);
                    string[] timeStr = new string[timePoints.Count];

                    for (int i = 0; i < timeStr.Length; i++)
                        timeStr[i] = ((float)(timePoints[i] / (1000.0 * 60.0))).ToString();//ms to min
                    //scanTimes.Add(time);//

                    foreach (CharacterizedPrecursor peptide in mixed.PeptideRatios.Keys)
                    {
                        scanTimes.Add(timeStr);
                        double[] arrayInt = new double[timePoints.Count];
                        for (int i = 0; i < arrayInt.Length; i++)
                            //arrayInt[i] = mixed.PeptideRatios[peptide].eCurvePerMs.InterpolateIntensity(timePoints[i]);
                            arrayInt[i] = mixed.PeptideRatios[peptide].eCurveCount.InterpolateIntensity(timePoints[i]);
                        intensities.Add(arrayInt);

                        ratios.Add(mixed.PeptideRatios[peptide].eCurvePerMs.Area);

                        names.Add(peptide.Peptide.Sequence + "  [" + Utilities.vsCSV.GetFileName_NoExtension(peptide.Sample.sSDF) + "]");
                    }//*/

                    StackPanel1.Children.Add(new SynthPeptideUC(intensities.ToArray(), scanTimes.ToArray(), names.ToArray(), Utilities.vsCSV.GetFileName_NoExtension(sample.sSDF) + " (Precursor " + mixed.MZ + ")"));
                    StackPanel1.Children.Add(new RatioUC(ratios.ToArray(), names.ToArray(), Utilities.vsCSV.GetFileName_NoExtension(sample.sSDF) + "  [" + mixed.MZ + "]"));
                }
            }

            //Export charts to result folder (only works if first render is done. Only called if the page is navigated to)
            //ExportToImageFile(PepIso.solver.OutputFolder + sample.Name + "_" + strKey + ".png");
            outputFile = PepIso.solver.OutputFolder + sample.Name + "_" + strKey + ".png";
        }

        private void ExportToImageFile(string outputFile)
        {
            Size size = new Size(StackPanel1.ActualWidth, StackPanel1.ActualHeight);
            RenderTargetBitmap bmpSrc = new RenderTargetBitmap((int)size.Width, (int)size.Height, 96, 96, PixelFormats.Pbgra32);

            DrawingVisual drawingvisual = new DrawingVisual();
            using (DrawingContext context = drawingvisual.RenderOpen())
            {
                context.DrawRectangle(new VisualBrush(StackPanel1), null, new Rect(new Point(), size));
                context.Close();
            }

            bmpSrc.Render(drawingvisual);

            PngBitmapEncoder encoder = new PngBitmapEncoder();
            encoder.Frames.Add(BitmapFrame.Create(bmpSrc));
                        
            StreamWriter sw = new StreamWriter(outputFile);
            encoder.Save(sw.BaseStream);
            sw.Close();
        }

        private void BtnSave_Click(object sender, RoutedEventArgs e)
        {
            ExportToImageFile(outputFile);
            txtSave.Text = outputFile;
            /*
            bool? dResult = false;
            Microsoft.Win32.SaveFileDialog save_dlg = new Microsoft.Win32.SaveFileDialog();
            save_dlg.Title = "Exporting chart as image";
            save_dlg.Filter = "Image file (*.png)|*.png";
            save_dlg.FileName = "image.png";
            //open_dlg.Multiselect = true;
            save_dlg.ValidateNames = true;
            save_dlg.OverwritePrompt = true;
            //save_dlg.
            dResult = save_dlg.ShowDialog();
            if (!string.IsNullOrEmpty(save_dlg.FileName))
                ExportToImageFile(save_dlg.FileName);//*/
        }
    }
}
