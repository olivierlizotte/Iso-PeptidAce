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

namespace PeptidAce.ModernUI.Content
{
    /// <summary>
    /// Interaction logic for About.xaml
    /// </summary>
    public partial class ResultsMixed : UserControl, IContent
    {
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
            Dictionary<string, string> parameters = GetParams(e.Source.OriginalString);
            if(parameters.ContainsKey("Sample") && parameters.ContainsKey("mzKey"))
            {
                foreach(Sample sample in PepIso.solver.mixedPrecursors.Keys)
                    if(sample.Name.CompareTo(parameters["Sample"]) == 0)
                        UpdateContent(sample, (double.Parse(parameters["mzKey"])));
            }
        }
        public void OnNavigatingFrom(NavigatingCancelEventArgs e)
        {
        }

        public ResultsMixed()
        {
            InitializeComponent();
        }

        public void UpdateContent(Sample sample, double mzKey)
        {
            foreach(MixedPrecursor mixed in PepIso.solver.mixedPrecursors[sample])
            {
                if (mixed.PeptideRatios.Count > 0 && mixed.MZ == mzKey)
                {
                    List<double> ratios = new List<double>();
                    List<string> names = new List<string>();

                    foreach (CharacterizedPrecursor peptide in mixed.PeptideRatios.Keys)
                    {
                        ratios.Add(mixed.PeptideRatios[peptide].eCurvePerMs.Area);
                        names.Add(peptide.Peptide.Sequence + "  [" + Utilities.vsCSV.GetFileName_NoExtension(peptide.Sample.sSDF) + "]");
                    }

                    StackPanel1.Children.Add(new RatioUC(ratios.ToArray(), names.ToArray(), Utilities.vsCSV.GetFileName_NoExtension(sample.sSDF) + "  [" + mixed.MZ + "]"));
                }
            }
        }
    }
}
