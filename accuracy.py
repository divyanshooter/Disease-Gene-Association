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

	ten_percent_mean_JI=0
	ten_percent_mean_SI=0
	ten_percent_mean_GI=0
	ten_percent_mean_CI=0
	for i in range(10):
		accuracy=[0,0,0,0]
		x = algorithm(Working_values[i],k) # DGA algorithm is run with working_values(modified dataset and choosen k)
		count_of_values = len(Missing_values[i])
		for j in x:
			if j[:2] in Missing_values[i]:
				for cnt in range(2,6):
					accuracy[cnt-2]+=j[cnt]
		for i in range(4):
			accuracy[i]=accuracy[i]/n   # mean of accuracy
		ten_percent_mean_JI += 0.1*accuracy[0] #contributing to 10% of accuracy as 10 crossover for Jaccard Index
		ten_percent_mean_SI += 0.1*accuracy[1] #contributing to 10% of accuracy as 10 crossover for Simpson Index
		ten_percent_mean_GI += 0.1*accuracy[2] #contributing to 10% of accuracy as 10 crossover for Geometric Index
		ten_percent_mean_CI += 0.1*accuracy[3] #contributing to 10% of accuracy as 10 crossover for Cosine Index
	return [ten_percent_mean_JI,ten_percent_mean_SI,ten_percent_mean_GI,ten_percent_mean_CI]
	
