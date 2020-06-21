# Importing necesary modules
import numpy as np
import pandas as pd
import networkx as nx
#import networkx.algorithms.approximation as nxaa
import matplotlib.pyplot as plt
import igraph as ig
import io, base64

def olke(seq, bases, keep_seq):
    # Defining weights of base-pairs for use in tables
    main_dict = {'GU':0, 'UG':0, 'AU':1, 'UA': 1, 'GC':2, 'CG':2}

    # Taking input
    input_arr = list(seq.upper().replace('T','U'))

    all_bases = bases + [i[::-1] for i in bases]
	# Storing into a dataframe
    df = pd.DataFrame(input_arr, columns=['base'])
    df.index.name = 'index'
	## print(df, end='\n\n\n')


	# Creating the table for storing weights of pairs acc. to main_dict
    length = len(input_arr)
	# table = -np.ones((length, length))

    # making the table
    table = np.array([[main_dict[input_arr[i] + input_arr[j]] 
					if input_arr[i]+input_arr[j] in all_bases else -1
					for j in range(length)]
					for i in range(length)])

    for i in range(length):
	    if keep_seq[i] == 'x' or keep_seq[i] == 'X':
		    table[i][:] = -1
		    table[:][i] = -1

    # storing in DFs and plotting the graph based on table
    store_array = []
    k = 0
    for i in range(length):
        for j in range(i, length):
          if table[i][j] != -1:
              store_array.append([i, j, table[i][j], k])
              k += 1

    edge_list = [tuple(i[:3]) for i in store_array]

    df2 = pd.DataFrame(store_array, columns=['from', 'to', 'weights', 'name'])

    graph = nx.Graph()
    graph.add_nodes_from(list(range(length)))
    graph.add_weighted_edges_from(edge_list, attr_dict={'weight': i[2] for i in store_array})

    # plt.figure(figsize=(12,8))
    # plt.title("The main graph with nodes corresponding to bases(starting with 0)")
    # nx.draw(graph, node_color="black", font_color="white", with_labels=True, pos=nx.circular_layout(graph))
    # nx.draw_networkx_edge_labels(graph, nx.circular_layout(graph), edge_labels={(i[0],i[1]) : i[3] for i in store_array}, font_color='red')
    # plt.show()




    # function for making a cycle graph
    from itertools import combinations

    def make_graph(df, length):
        inters_list = pd.DataFrame(columns=['from', 'to'])
        edge_combs = list(combinations(df['name'], 2))
        # print('Length', len(edge_combs))
        count = 0
        for i in edge_combs:
          x = df.loc[df['name'] == i[0]]
          y = df.loc[df['name'] == i[1]]

          if (x['from'].item() < y['from'].item()) and (y['from'].item() < x['to'].item()) and (x['to'].item() < y['to'].item()):
              inters_list.loc[count] = {'from': x['name'].item(), 'to': y['name'].item()}
              count += 1
                
          elif (y['from'].item() < x['from'].item()) and (x['from'].item() < y['to'].item()) and (y['from'].item() < x['to'].item()):
              inters_list.loc[count] = {'from': y['name'].item(), 'to': x['name'].item()}
              count += 1

          elif (x['to'].item() == y['to'].item()):
              inters_list.loc[count] = {'from': x['name'].item(), 'to': y['name'].item()}
              count += 1

          elif (x['from'].item() == y['to'].item()):
              inters_list.loc[count] = {'from': x['name'].item(), 'to': y['name'].item()}
              count += 1

          elif (x['to'].item() == y['from'].item()):
              inters_list.loc[count] = {'from': x['name'].item(), 'to': y['name'].item()}
              count += 1

          elif (x['from'].item() == y['from'].item()):
              inters_list.loc[count] = {'from': x['name'].item(), 'to': y['name'].item()}
              count += 1
            

        return inters_list



    # The table for identifying the presence of stacks in a given sequence
    def make_stack_table(df):
        table = -np.ones((length, length))
        for i in range(len(df)):
          table[df['from'][i]][df['to'][i]] = 1

        for x in range(length-1):
          for y in range(x+1, length):
              if table[x][y] == 1:
                  x_s = x + 1
                  for z in range(y-1,0, -1):
                      if table[x_s][z] == 1:
                          table[x][y] += 1
                          x_s += 1
                      else:
                          break

        return table


    # Counting the number of stacks after identifying them
    def count_stacks(stack_table):
        count = 0
        for x in range(length-1):
          for y in range(x+1, length):
              if stack_table[x][y] > 1:
                  count += 1

        return count


    # Constructing an array of stacks for a given config
    def stacks(st):
        stack_array = []
        for x in range(length-1):
          for y in range(x+1, length):
              if st[x][y] > 1:
                  if x==0 or y==(length-1):
                      stack_array.append(st[x][y]) 
                  elif st[x-1][y+1] < st[x][y]:
                      stack_array.append(st[x][y])
                        
        return sorted(stack_array,reverse=True)


    # Compare two possible secondary structures to find which one has more consecutive stacks
    def compare_stack_array(st1, st2):
        length = min(len(st1), len(st2))
        for i in range(length):
          if st1[i] > st2[i]:
              return 0
          elif st1[i] < st2[i]:
              return 1
        return 2


    # If the number of consecutive stacks turns out to the same, we compare the total base-pair energies of both and choose 
    # whichever has less Gibbs energy(based on energy_dict)

    energy_dict = {'GU':24.5, 'UG':24.5, 'AU':39.0, 'UA': 39.0, 'GC':78.8, 'CG':78.8}
    def compare_energy(st1, st2):
        if st1 < st2:
          return 0
        elif st1 > st2:
          return 1
        else:
          return 2


    # Stack table, each value corresponding to stacking energy upon one base pair stacking over the other
    stacking_energy_table = np.array([[-0.9, -1.8, -2.3, -1.1, -1.1, -0.8],
                                    [-1.7, -2.9, -3.4, -2.3, -2.1, -1.4],
                                    [-2.1, -2.0, -2.9, -1.8, -1.9, -1.2],
                                    [-0.9, -1.7, -2.1, -0.9, -1.0, -0.5],
                                    [-0.5, -1.2, -1.4, -0.8, -0.4, -0.2],
                                    [-1.0, -1.9, -2.1, -1.1, -1.5, -0.4]])
     

    # To interpret the stack value in the table
    def find_row_or_col(nxt):
        if nxt == 'AU':
          return 0
        elif nxt == 'CG':
          return 1
        elif nxt == 'GC':
          return 2
        elif nxt == 'UA':
          return 3
        elif nxt == 'GU':
          return 4
        elif nxt == 'UG':
          return 5


    # Calculating total stack energy for a given structure 
    def cal_stack_energy(table):
        total_energy = 0
        for p in range(length-1):
          for q in range(p+1, length):
              if table[p][q] > 1:
                  ## print("p and q", p, q)
                  ind1 = find_row_or_col(input_arr[p] + input_arr[q])
                  ind2 = find_row_or_col(input_arr[p+1] + input_arr[q-1])
                  total_energy += stacking_energy_table[ind1][ind2]
        return total_energy


    # Calculating the bond energy for all base-pairs
    def cal_bond_energy(table):
        total_energy = 0
        for p in range(length-1):
          for q in range(p+1, length):
              if table[p][q] > 0:
                  ## print("p and q", p, q)
                  ind = input_arr[p] + input_arr[q]
                  total_energy += energy_dict[ind]
        return total_energy



    def cal_olke_energy(base_pairs):
        total_energy = main_dict[input_arr[base_pairs[0][0]] + input_arr[base_pairs[0][1]]]
        for i in range(len(base_pairs)-1):
          down_start = base_pairs[i][0] + 1
          up_start = base_pairs[i+1][0]

          down_end = base_pairs[i][1] - 1
          up_end = base_pairs[i+1][1]

          total_energy += main_dict[input_arr[up_start] + input_arr[up_end]]
            

          if up_start == down_start:
              # bulge towards right
              if down_end != up_end:

                  if (down_end - up_end) <= 2:
                      total_energy -= 2
                  elif (down_end - up_end) == 3:
                      total_energy -= 3
                  elif (down_end - up_end) <= 15:
                      total_energy -= 5
                  else:
                      total_energy -= 6
                    
          else:
              # bulge towards left
              if down_end == up_end:

                  if (up_start - down_start) <= 2:
                      total_energy -= 2
                  elif (up_start - down_start) == 3:
                      total_energy -= 3
                  elif (up_start - down_start) <= 15:
                      total_energy -= 5
                  else:
                      total_energy -= 6
                
              elif down_end != up_end:
                  # interior loops
                  if up_start < down_end: 
                      if (down_end - up_end) + (up_start - down_start) == 2:
                          total_energy -= 4
                      elif (down_end - up_end) + (up_start - down_start) == 3:
                          total_energy -= 5
                      elif (down_end - up_end) + (up_start - down_start) <= 7:
                          total_energy -= 6
                      else:
                          total_energy -= 7

                  # hairpin loops
                  else:
                      hair_pin_size = (down_end - down_start) - 1
                      if hair_pin_size == 3:
                          total_energy -= 5
                      elif hair_pin_size <= 7:
                          total_energy -= 6
                      else:
                          total_energy -= 7

          if i == len(base_pairs)-2:
              hair_pin_size = (down_end - down_start) - 1
              # print(down_end, down_start)
              if hair_pin_size == 3:
                  total_energy -= 5
              elif hair_pin_size <= 7:
                  total_energy -= 6
              else:
                  total_energy -= 7      
                      
        return total_energy


        


    # making the circle graph using the make_graph function
    cir_graph = nx.Graph()
    cir_graph.add_nodes_from(list(range(k)))
    df2 = df2.loc[((df2['to'] - df2['from']) > 3)]
    cir_edge_frame = make_graph(df2, length)
    cir_edge_list = list(zip(cir_edge_frame['from'], cir_edge_frame['to']))
    cir_graph.add_edges_from(cir_edge_list)


    # Converting the networkx graph to a python-igraph due to its better support for finding all maximum independent vertex sets
    # for a graph
    g = ig.Graph(directed=False)
    g.add_vertices(cir_graph.nodes())
    g.add_edges(cir_graph.edges())
    mis_max = ig.Graph.largest_independent_vertex_sets(g) 
    # print("Total maximal independent sets:", len(mis_max))


    # Imposing restrictions on the MISs one after the other.
    # 1) We choose the structures with the maximum number of stacking regions
    # 2) If there are several possible such structures, we then check which structures have the most no of consecutive stacks
    # 3) If there are several possible such structures, we then compare the energies of the stacks based on the 'stacking energy table'
    # 4) If there are several possible such structures, we then compare the individual bond energies.

    # After these operations, we plot the structures obtained. These we say are the best possible structures. 
    max_stacks = 0
    best_table = []
    best_table_df = []
    best_stack_array = [0]
    best_stack_energy = best_bond_energy = best_olke_energy = 0
    k = 1
    for i in mis_max:
        req_df = df2[df2['name'].isin(i)].reset_index()
        # con_graph = nx.Graph()
        # con_graph.add_nodes_from(list(range(length)))
        # con_graph.add_edges_from(list(zip(req_df['from'], req_df['to'])))

        # nx.draw(con_graph, node_color="black", font_color="white", with_labels=True, pos=nx.circular_layout(con_graph))
        # con_edge = dict(zip(zip(req_df['from'], req_df['to']), req_df['name']))
        # nx.draw_networkx_edge_labels(con_graph, nx.circular_layout(con_graph), edge_labels=con_edge, font_color='red')


        stack_table = make_stack_table(req_df)
        total_stacks = count_stacks(stack_table)
        stack_array = stacks(stack_table)
        stack_energy = cal_stack_energy(stack_table)
        bond_energy = cal_bond_energy(stack_table)
        okle_energy = cal_olke_energy(sorted(list(zip(req_df['from'], req_df['to']))))

        # Filtering out sets with the maximum stacking regions
        if total_stacks == max_stacks:
          array_flag = compare_stack_array(stack_array, best_stack_array)
          if array_flag == 0:
              # If more stacks in a stacking region, discard previous and keep this
              best_table_df.clear()
              best_table_df.append(req_df)
              best_stack_array = stack_array
          elif array_flag == 2:
              # compare their olke energies
              olke_flag = compare_energy(okle_energy, best_olke_energy)
              if olke_flag == 1:
                  best_table_df.clear()
                  best_table_df.append(req_df)
                  best_olke_energy = okle_energy

              elif olke_flag == 2:
                  # If equivalent check for stack-energy
                energy_flag = compare_energy(stack_energy, best_stack_energy)
                if energy_flag == 0:
                  best_table_df.clear()
                  best_table_df.append(req_df)
                  best_stack_energy = stack_energy

                



        if total_stacks > max_stacks:
          best_table_df.clear(); best_stack_array.clear()
          best_table_df.append(req_df)
          best_stack_array, max_stacks, best_stack_energy, best_bond_energy = stack_array, total_stacks, stack_energy, bond_energy
          best_olke_energy = okle_energy
             
        k += 1
        #plt.show()

    # print("Olke", best_olke_energy)
    # print("Total sets after all restrictions:", len(best_table_df))
    # print(best_stack_energy)

    def fig_to_base64(fig):
        img = io.BytesIO()
        fig.savefig(img, format='png', bbox_inches='tight')
        img.seek(0)
        return base64.b64encode(img.getvalue())
    
    output_graphs = []
    for req_df in best_table_df:
        fig = plt.figure(figsize=(4,3))
        con_graph = nx.Graph()
        con_graph.add_nodes_from(list(range(length)))
        con_graph.add_edges_from(list(zip(req_df['from'], req_df['to'])))

        nx.draw(con_graph, node_color="black", font_color="white", with_labels=True, pos=nx.circular_layout(con_graph))
        con_edge = dict(zip(zip(req_df['from'], req_df['to']), req_df['name']))
        nx.draw_networkx_edge_labels(con_graph, nx.circular_layout(con_graph), edge_labels=con_edge, font_color='red')
        encoded = fig_to_base64(fig)
        # my_html = '<img src="data:image/png;base64, {}">'.format(encoded.decode('utf-8'))
        my_html = encoded.decode('utf-8')
        output_graphs.append(my_html)
    



    # Expressing the results in a dot-bracket notation
    dot_bracket_array = []
    for req_df in best_table_df:
        good_df = req_df[['from', 'to']]
        final_edges = list(zip(good_df['from'], good_df['to']))

        dot_bracket = ['.' for _ in range(length)]
        for i in final_edges:
          dot_bracket[i[0]], dot_bracket[i[1]] = '(', ')'
          
        dot_bracket_array.append(''.join(dot_bracket))

    return dot_bracket_array, output_graphs


