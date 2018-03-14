# filter rows in pandas:
# sum  all ground-state-counts (#ground) for run 1:
np.sum(df[(df['run'] == 1)]['#ground'])
