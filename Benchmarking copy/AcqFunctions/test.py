"""
Created on Wed Aug 20 14:06:50 2025

@author: josephbailey
"""

import pandas as pd
import numpy as np
def process_bayesopt_data(csv_file):
    # Read the CSV file
    df = pd.read_csv(csv_file)
    
    # Initialize matrices - we have runs 0-4 (5 runs total)
    num_runs = 5
    num_trials = 20
    
    # Create matrices for binding and folding scores (direct min values)
    binding_matrix = [[0.0 for _ in range(num_trials)] for _ in range(num_runs)]
    folding_matrix = [[0.0 for _ in range(num_trials)] for _ in range(num_runs)]
    
    # Process each run (0 to 4)
    for run in range(num_runs):
        # Filter data for current run (0-4)
        run_data = df[df['Run'] == run]
        print(f"Run {run}: Found {len(run_data)} data points")
        
        # Create a dictionary to store trial data for easy lookup
        trial_dict = {}
        for _, row in run_data.iterrows():
            trial = row['Points Evaluated']
            if 1 <= trial <= 20:  # Only consider trials 1-20
                trial_dict[trial] = {
                    'binding': row['Binding ddG'],
                    'folding': row['Folding ddG']
                }
                print(f"  Trial {trial}: Binding={row['Binding ddG']:.6f}, Folding={row['Folding ddG']:.6f}")
        
        # Initialize current minimums with first data point from the entire dataset
        current_min_binding = float(df['Binding ddG'].iloc[0])
        current_min_folding = float(df['Folding ddG'].iloc[0])
        print(f"Run {run}: Initial minimums - Binding: {current_min_binding}, Folding: {current_min_folding}")
        
        # Process each trial in order (1 to 20) - NOTE: Fixed indentation
        for trial_num in range(1, num_trials + 1):
            if trial_num in trial_dict:
                # Trial exists in data - get actual values
                current_binding = float(str(trial_dict[trial_num]['binding']).strip())
                current_folding = float(str(trial_dict[trial_num]['folding']).strip())
                
                # Update current minimums if we found a better value
                if current_binding < current_min_binding:
                    current_min_binding = current_binding
                if current_folding < current_min_folding:
                    current_min_folding = current_folding
                
                print(f"  Trial {trial_num}: Updated mins - Binding: {current_min_binding:.6f}, Folding: {current_min_folding:.6f}")
            else:
                print(f"  Trial {trial_num}: No data, keeping previous mins - Binding: {current_min_binding:.6f}, Folding: {current_min_folding:.6f}")
            
            # Store the current minimum (up to this trial)
            binding_matrix[run][trial_num-1] = current_min_binding
            folding_matrix[run][trial_num-1] = current_min_folding
    
    # Convert to DataFrames for better presentation
    binding_df = pd.DataFrame(binding_matrix, 
                             index=[f'Run {i}' for i in range(num_runs)],
                             columns=[f'Trial {i+1}' for i in range(num_trials)])
    
    folding_df = pd.DataFrame(folding_matrix,
                             index=[f'Run {i}' for i in range(num_runs)],
                             columns=[f'Trial {i+1}' for i in range(num_trials)])
    
    return binding_df, folding_df




# Example usage
if __name__ == "__main__":
    binding_matrix, folding_matrix = process_bayesopt_data('/Users/nathanphan/Desktop/Projects/ncPPI/Benchmarking/AcqFunctions/08_24_u_UCB.csv')
    
    # binding_matrix = round_min(binding_matrix)
    # folding_matrix = round_min(folding_matrix)
    print("\nBinding ddG Matrix (All Minimum Values):")
    print(binding_matrix)
    print("\n" + "="*80 + "\n")
    print("Folding ddG Matrix (All Minimum Values):")
    print(folding_matrix)
    
    # Optional: Save to CSV files
    binding_matrix.to_csv('U_UCB_Binding.csv')
    folding_matrix.to_csv('U_UCB_Folding.csv')