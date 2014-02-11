using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace PeptidAce.ModernUI.Content
{
    /// <summary>
    /// Interaction logic for SettingsAppearance.xaml
    /// </summary>
    public partial class SettingsAppearance : UserControl
    {
        public SettingsAppearance()
        {
            InitializeComponent();

            // create and assign the appearance view model
            //this.DataContext = new SettingsAppearanceViewModel();
            PrecMassTol.Text = ModernUI.Content.PepIso.precursorMassTolPPm.ToString();
            ProdMassTol.Text = ModernUI.Content.PepIso.productMassTolPPm.ToString();
            minNbFragment.Text = ModernUI.Content.PepIso.nbMinFragments.ToString();
            maxNbFragment.Text = ModernUI.Content.PepIso.nbMaxFragments.ToString();
        }

        private void TextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            try
            {
                ModernUI.Content.PepIso.precursorMassTolPPm = double.Parse(PrecMassTol.Text);
            }
            catch(Exception)
            {
            }            
        }

        private void TextBox_TextChanged_1(object sender, TextChangedEventArgs e)
        {
            try
            {
                ModernUI.Content.PepIso.productMassTolPPm = double.Parse(ProdMassTol.Text);
            }
            catch (Exception)
            {
            }            
        }

        private void TextBox_TextChanged_2(object sender, TextChangedEventArgs e)
        {
            try
            {
                ModernUI.Content.PepIso.nbMinFragments = int.Parse(minNbFragment.Text);
            }
            catch (Exception)
            {
            }
        }
        private void TextBox_TextChanged_3(object sender, TextChangedEventArgs e)
        {
            try
            {
                ModernUI.Content.PepIso.nbMaxFragments = int.Parse(maxNbFragment.Text);
            }
            catch (Exception)
            {
            }
        }
    }
}
