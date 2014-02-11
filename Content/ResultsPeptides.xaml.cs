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
    public partial class ResultsPeptides : UserControl, IContent
    {
        private string outputFile = null;
        public static Dictionary<string, string> GetParams(string uri)
        {
            string[] splits = uri.Split(new char[] { '?', '&' });
            Dictionary<string, string> parameters = new Dictionary<string, string>();
            foreach (string str in splits)
            {
                string[] pairs = str.Split('=');
                if (pairs.Length > 1)
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
                foreach (double key in PepIso.solver.characterizedPeptides.Keys)
                    if (key.ToString().CompareTo(parameters["mzKey"]) == 0)
                        UpdateContent(parameters["Sample"], key);
            }
        }
        public void OnNavigatingFrom(NavigatingCancelEventArgs e)
        {
        }

        public ResultsPeptides()
        {
            InitializeComponent();
        }

        public void UpdateContent(string sampleName, double key)
        {
            //foreach(double key in PepIso.solver.characterizedPeptides.Keys)
            //{
            foreach (Sample sample in PepIso.solver.characterizedPeptides[key].Keys)
            {
                if (sample.Name.CompareTo(sampleName) == 0)
                {
                    CharacterizedPrecursor cPrec = PepIso.solver.characterizedPeptides[key][sample];

                    List<string[]> scanTimes = new List<string[]>();
                    List<double> timeArray = Utilities.ElutionCurve.GetTimePoints(64, true, cPrec.eCurve.interpolatedTime, cPrec.eCurve.interpolatedIntensityCount);
                    List<string> time = new List<string>(timeArray.Count);
                    foreach (double timeP in timeArray)
                        time.Add(((float)(timeP / (1000.0 * 60.0))).ToString());
                    scanTimes.Add(time.ToArray());

                    List<double[]> intensities = new List<double[]>();
                    double[] arrayInt = new double[timeArray.Count];
                    for (int i = 0; i < arrayInt.Length; i++)
                        arrayInt[i] = cPrec.eCurve.InterpolateIntensity(timeArray[i]);
                    intensities.Add(arrayInt);
                    //intensities.Add(cPrec.eCurve.intensityCount.ToArray());
                    
                    string[] names = new string[1];
                    names[0] = cPrec.Peptide.Sequence + "  [" + Utilities.vsCSV.GetFileName_NoExtension(cPrec.Sample.sSDF) + "]";
                    StackPanel1.Children.Add(new SynthPeptideUC(intensities.ToArray(), scanTimes.ToArray(), names, Utilities.vsCSV.GetFileName_NoExtension(sample.sSDF) + " (Precursor " + key + ")"));
                }
            }
            //}
            //Export charts to result folder (only works if first render is done. Only called if the page is navigated to)
            //ExportToImageFile(PepIso.solver.OutputFolder + sample.Name + "_" + strKey + ".png");
            outputFile = PepIso.solver.OutputFolder + sampleName + "_" + key.ToString() + ".png";
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
