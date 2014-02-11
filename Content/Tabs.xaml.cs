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
using FirstFloor.ModernUI.Windows;
using FirstFloor.ModernUI.Windows.Navigation;
using FirstFloor.ModernUI.Presentation;
using PeptidAce.Iso.Structures;
using PeptidAce.Iso.Methods;

namespace PeptidAce.ModernUI.Content
{
    /// <summary>
    /// Interaction logic for Tabs.xaml
    /// </summary>
    public partial class Tabs : UserControl, IContent
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

        //public bool isUpToDate = false;

        public Iso.Methods.PositionnalIsomerSolver solverPTR;
        public void OnNavigatedTo(NavigationEventArgs e)
        {
            //if (!isUpToDate)
            if (solverPTR != PepIso.solver)
            {
                solverPTR = PepIso.solver;
                //isUpToDate = true;
                ModernTabList.Links.Clear();
                ModernTabList.SelectedSource = null;
                Dictionary<string, string> parameters = GetParams(e.Source.OriginalString);
                if (parameters.ContainsKey("Sample"))
                {
                    string sampleName = parameters["Sample"];
                    foreach (double mzKey in solverPTR.characterizedPeptides.Keys)
                    {
                        if (solverPTR.characterizedPeptides[mzKey].Count > 0)
                        {
                            //Mixed Deconvoluted curves
                            foreach (Sample sample in solverPTR.mixedPrecursors.Keys)
                            {
                                string name = sample.Name;
                                if (name.CompareTo(sampleName) == 0)
                                {
                                    bool found = false;
                                    foreach (MixedPrecursor mp in solverPTR.mixedPrecursors[sample])
                                        if (mp.MZ == mzKey && mp.PeptideRatios.Count > 0)
                                            found = true;
                                    if (found)
                                    {
                                        Uri source = new Uri("/Content/ResultsDeconvoluted.xaml?mzKey=" + mzKey.ToString() + "&Sample=" + name, UriKind.Relative);
                                        ModernTabList.Links.Add(new Link
                                        {
                                            DisplayName = ((float)mzKey).ToString(),
                                            Source = source
                                        });
                                        if (ModernTabList.SelectedSource == null)
                                            ModernTabList.SelectedSource = source;
                                    }
                                    break;
                                }
                            }

                            //Mixed Results
                            /*foreach (Sample sample in PepIso.solver.mixedPrecursors.Keys)
                            {
                                string name = sample.Name;
                                ModernTabList.Links.Add(new Link
                                {
                                    DisplayName = name + " Curves",
                                    Source = new Uri("/Content/ResultsMixed.xaml?mzKey=" + mzKey + "&Sample=" + name, UriKind.Relative)
                                });
                            }//*/
                        }
                    }
                }
            }
        }
        public void OnNavigatingFrom(NavigatingCancelEventArgs e)
        {
        }

        public Tabs()
        {
            InitializeComponent();
        }
    }
}
