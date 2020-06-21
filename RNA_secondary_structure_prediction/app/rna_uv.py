# Importing necesary modules
import numpy as np
import pandas as pd
import networkx as nx
#import networkx.algorithms.approximation as nxaa
import matplotlib.pyplot as plt
import igraph as ig
import time
import io, base64


def u_v(seq, bases, keep_seq):
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

	# Deleting the table as we don't need it any further	
	del(table)

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
		#maxm = int(input("Enter the maximum value: "))
		#minm = int(input("Enter the minimum value: "))
		A = 1; B = 1
		h = lambda x: 1 if x==0 else 0

		n = len(df)
		# # print('len', n)
		U = [-5 for i in range(n)]
		V = [0 for i in range(n)]
		

		for time_step in range(800):
			# print("Time step:", time_step)
			stop_flag = 0
			# # print("iter", time_step)
			# V = [1 if U[i] > 0 else 0 for i in range(n)]
			for i in range(n):
				V[i] = 1 if U[i]>0 else 0
			for i in range(n):
				x1 = df.loc[i]['from']
				y1 = df.loc[i]['to']
				## print(x1, y1)	

				# V = [1 if U[i] > 0 else 0 for i in range(n)]
				delta_U_i = 0

				term = 0
				distance_i = min(abs(x1-y1), abs(length+x1-y1))
				# # print("distance_i", distance_i)
				for j in range(n):
					if j==i:
						continue
					d = 0
					x2 = df.loc[j]['from']
					y2 = df.loc[j]['to']

					if ((x1 < x2) and (x2 < y1) and (y1 < y2)) or ((x2 < x1) and (x1 < y2) and (y2 < y1)) or ((y1 == y2) or (x1 == y2) or (y1 == x2) or (x1 == x2)):
						d = 1	
					
					# V[j] = 1 if U[j] > 0 else 0

					term += (d*(1-V[j])) 
					# # print("j is {} and term is {}".format(j, term))
					
				delta_U_i = A*term*(1-V[i])/distance_i - B*h(term)*V[i]
				# # print('i is {} and delta_U_i is {}'.format(df.loc[i]['name'], delta_U_i))

				# # print(delta_U_i)
				if delta_U_i != 0:
					stop_flag = 1

				U[i] += delta_U_i
				# V[i] = 1 if U[i] > 0 else 0
			# print("U is", U)
			# print("V is", V)

			if stop_flag == 0:
				# print("Stopped in iter_no :", time_step+1)
				# print('\n\n')
				break


		# print("V is",V)
		# print('\n\n\n')
		edges = [not bool(i) for i in V]
		edges = pd.Series(edges)	
		# # print(df.loc[edges, :])
		return df.loc[edges, :]
		



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
	          # # print(down_end, down_start)
	          if hair_pin_size == 3:
	              total_energy -= 5
	          elif hair_pin_size <= 7:
	              total_energy -= 6
	          else:
	              total_energy -= 7      
	                  
	    return total_energy



	# making the circle graph using the make_graph function
	cir_graph = nx.Graph()
	cir_graph.add_nodes_from(list(range(length)))
	df2 = df2.loc[((df2['to'] - df2['from']) > 3)]
	df2 = df2.reset_index(drop=True)
	cir_edge_frame = make_graph(df2, length)
	cir_edge_frame = cir_edge_frame.reset_index(drop=True)
	# print(cir_edge_frame)
	mis_max = []
	mis_max.append(list(cir_edge_frame['name']))
	cir_edge_list = list(zip(cir_edge_frame['from'], cir_edge_frame['to']))
	cir_graph.add_edges_from(cir_edge_list)



	def fig_to_base64(fig):
		img = io.BytesIO()
		fig.savefig(img, format='png', bbox_inches='tight')
		img.seek(0)
		return base64.b64encode(img.getvalue())

	output_graphs = []
	fig = plt.figure(figsize=(4,3))
	# plt.title("The graph after applying U,V conditions")
	nx.draw(cir_graph, node_color="black", font_color="white", with_labels=True, pos=nx.circular_layout(cir_graph))
	c = cir_edge_frame
	nx.draw_networkx_edge_labels(cir_graph, nx.circular_layout(cir_graph), 
			edge_labels= {(c.loc[i]['from'], c.loc[i]['to']) : c.loc[i]['name'] for i in range(len(c))}, font_color='red')
	encoded = fig_to_base64(fig)
    # my_html = '<img src="data:image/png;base64, {}">'.format(encoded.decode('utf-8'))
	my_html = encoded.decode('utf-8')
	output_graphs.append(my_html)



	# Expressing the results in a dot-bracket notation
	dot_bracket_array = []
	good_df = cir_edge_frame[['from', 'to']]
	final_edges = list(zip(good_df['from'], good_df['to']))
	dot_bracket = ['.' for _ in range(length)]
	for i in final_edges:
		dot_bracket[i[0]], dot_bracket[i[1]] = '(', ')'

	dot_bracket_array.append(''.join(dot_bracket))

	return dot_bracket_array, output_graphs

	# print("\n\nDot bracket notation:", ''.join(dot_bracket))
	# print("\n\nTotal stability number: ", cal_olke_energy(cir_edge_list))
	# print("\n\nTotal time taken: {} seconds".format(time.time() - start))