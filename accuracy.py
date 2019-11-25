def ten_crossover(values,algorithm,k):
#values of dataset, DGA algorithm, k - number of nearest neighbour

	Missing_values=[[] for i in range(10)]
	Working_values=[[] for i in range(10)]
	for i in values:
		x=random.randint(0,9)
		for j in range(10):
			if j==x:
				Missing_values[j].append(i)
			else:
				Working_values[j].append(i)

	ten_percent_mean=0
	for i in range(10):
		accuracy=0
		x = algorithm(Working_values,k) # DGA algorithm is run with working_values(modified dataset and choosen k)
		count_of_values = len(Missing_values)
		for j in x:
			if j[:2] in Missing_values:
				accuracy+=j[2]
		accuracy=accuracy/n # mean of accuracy
		ten_percent_mean += 0.1*accuracy #contributing to 10% of accuracy as 10 crossover

	return ten_percent_mean