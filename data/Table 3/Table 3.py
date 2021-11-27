from VGsim import Simulator
import time
import os 

def get_data(num_pop, total_migration_rate):
    simulator = Simulator(number_of_sites=2, populations_number=num_pop, number_of_susceptible_groups=3, seed=2023)
    simulator.set_transmission_rate(2.5)
    simulator.set_transmission_rate(4.0, haplotype='GG')
    simulator.set_recovery_rate(0.9)
    simulator.set_sampling_rate(0.1)
    simulator.set_migration_probability(total_probability=total_migration_rate)
    simulator.set_population_size(int(2*10**9/num_pop))
    simulator.set_susceptibility(1.0, susceptibility_type=0)
    simulator.set_susceptibility(0.0, susceptibility_type=1)
    simulator.set_susceptibility(0.5, susceptibility_type=1, haplotype='C*')
    simulator.set_susceptibility(1.0, susceptibility_type=1, haplotype='G*')
    simulator.set_susceptibility(0.0, susceptibility_type=2)
    simulator.set_susceptibility_type(1)
    simulator.set_susceptibility_type(2, haplotype='C*')
    simulator.set_susceptibility_type(2, haplotype='G*')
    simulator.set_immunity_transition(0.01, target=0)
    simulator.set_lockdown([0.1, 0.02, 0.01])
    time_1 = time.time()
    simulator.simulate(iterations=100000000)
    time_2 = time.time()
    
    return [str(round(time_2 - time_1, 1)), str(round(simulator.get_proportion()*100, 2))]

population_numbers = [2, 5, 10, 20, 50, 100]
total_migration_rates = [0.001, 0.002, 0.005, 0.01, 0.1]

with open('table.txt', 'w') as table:
    table.write('\\begin{table}[ht]\n')
    table.write('\t\\resizebox{\\linewidth}{!}{\n')
    table.write('\t\t\\begin{tabular}{|c||c|c|c|c|c|c|}\n')
    table.write('\t\t\t\\hline\n')
    table.write('\t\t\t\t\\multirow{2}{*}{\\shortstack[c]{Cumulative\\\\migration\\\\probability $M$}} & \\multicolumn{6}{c|}{\\makecell{Number of demes $K$\\\\\\,}} \\\\\n')
    table.write('\t\t\t\\cline{2-7}\n')
    table.write('\t\t\t\t& 2 & 5 & 10 & 20 & 50 & 100\\\\\n')
    table.write('\t\t\t\\hline\n')
    for rate in total_migration_rates:
        c = []
        for number in population_numbers:
            c.append(get_data(number, rate))
        table.write('\t\t\t\t'+str(rate)+' & '+c[0][0]+'s & '+c[1][0]+'s & '+c[2][0]+'s & '+c[3][0]+'s & '+c[4][0]+'s & '+c[5][0]+'s \\\\')
        table.write('\t\t\t\t& '+c[0][1]+'\\% & '+c[1][1]+'\\% & '+c[2][1]+'\\% & '+c[3][1]+'\\% & '+c[4][1]+'\\% & '+c[5][1]+'\\% \\\\\n')
        table.write('\t\t\t\\hline\n')
    table.write('\t\t\\end{tabular}\n')
    table.write('\t}\n')
    table.write('\\end{rable}\n')
