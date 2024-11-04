import plotly.express as px
import pandas as pd
import sys
from math import pi 

def main():
    if len(sys.argv) != 2:
        print("Usage: python visualization.py <prediction_file>")
        sys.exit(1)

    data_csv = pd.read_csv(sys.argv[1], sep='\s+')
    
    # Keep only the rows for which the column theta_i is equal to 0
    data_csv = data_csv[abs(data_csv['theta_i']-pi) < 1e-2]

    # The column id_man_comb is the prediction of the type of maneouvre and depends 
    # on the values of columns theta_f, alpha_m and alpha_f. Please draw a 3d plot

    fig1 = px.scatter_3d(data_csv, x='theta_f', y='alpha_m', z='alpha_f', color='id_man_comb',
                         labels={'theta_f': 'Theta_f', 'alpha_m': 'Alpha_m', 'alpha_f': 'Alpha_f'})
    fig1.update_layout(coloraxis_colorbar=dict(title='Manoeuvre Type'))
    fig1.show()

    # Project the values into a 3D space where the z is the id_man_comb and the x and y are 
    # theta_if and alpha_f respectively 

    fig2 = px.scatter_3d(data_csv, x='id_man_comb', y='alpha_f', z='theta_f',
                         labels={'id_man_comb': 'Manoeuvre Type', 'alpha_f': 'Alpha_f', 'theta_f': 'Theta_i'})
    fig2.show()

if __name__ == "__main__":
    main()