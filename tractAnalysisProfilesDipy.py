#!/usr/bin/env python3

def generateProfilesDatatype(tract_data,out_path):

	import pandas as pd

	# generate pandas dataframe
	out_df = pd.DataFrame(tract_data)

	# save dataframe
	out_df.to_csv(out_path,index=False)

def generateImagesDatatype(tract_json,tract_data,tract_name,measure_name,out_path):

	import json
	import numpy as np
	import seaborn as sns
	import matplotlib.pyplot as plt

	# if measure = odi, ndi, isovf, fa, md, ad, rd, ga; scale = 0 1, ticks = [0 0.25 0.5 0.75]; else scale = 0 2, ticks = [0 0.5 1 1.5]
	if measure_name in ['fa','ga','ndi','odi','isovf']:
		y_scale = [0,1]
		y_ticks = [0,0.25,0.5,0.75]
		dm_tag = "(unitless)"
	else:
		y_scale = [0,2]
		y_ticks = [0,0.5,1,1.5]
		dm_tag = "(um^2/msec)"

	x_ticks = [0,len(tract_data)]
	xticks_labels = ['Tract RAS','Tract LPI']

	x_label = 'Location on Tract'
	y_label = measure_name + '\n'+dm_tag

	fig, ax = plt.subplots(1)
	ax.plot(tract_data)
	ax.set_ylabel(y_label)
	ax.set_xlabel(x_label)
	ax.set_ylim(y_scale)
	ax.set_yticks(y_ticks)
	ax.set_xticks(x_ticks)
	ax.set_xticklabels(xticks_labels)
	ax.set_title(tract_name)

	sns.despine(fig)

	fig.savefig(out_path)

	plt.close()

	tmp = {}
	tmp['filename'] = 'images/'+out_path.split('/')[-1]
	tmp['name'] = tract_name
	tmp['desc'] = measure_name
	tract_json['images'] = np.append(tract_json['images'],tmp)

	return tract_json

def generateTractmeasuresDatatype(subjectID,names,measure_names,tract_data,out_path):

	import pandas as pd

	# create dataframe
	out_df = pd.DataFrame([])

	# append data to output df
	for bundles in names:
		tmp_df = pd.DataFrame(tract_data[bundles])
		tmp_df['subjectID'] = [ subjectID for f in range(len(tmp_df)) ]
		tmp_df['structureID'] = [ bundles for f in range(len(tmp_df)) ]
		tmp_df['nodeID'] = [ f+1 for f in range(len(tmp_df)) ]
		tmp_df = tmp_df[['subjectID','structureID','nodeID'] + [ f+'_mean' for f in measure_names ] ]
		tmp_df.columns = tmp_df.columns.str.rstrip('_mean')
		out_df = pd.concat([out_df,tmp_df])

	# output
	out_df.to_csv(out_path,index=False)

def computeTractProfiles(subjectID,reference_anat_path,streamlines_path,classification_path,measures_path,n_points,out_path):

	import os, sys
	import numpy as np, nibabel as nib, scipy.io as sio, pandas as pd
	import dipy.stats.analysis as dsa
	import dipy.tracking.streamline as dts
	from dipy.segment.clustering import QuickBundles
	from dipy.segment.metric import (AveragePointwiseEuclideanMetric,
	                                 ResampleFeature)
	from dipy.io.streamline import load_tractogram

	# load reference anatomy (dwi)
	print('loading reference anatomy')
	ref_anat = nib.load(reference_anat_path)

	# load tractogram
	print('loading tractogram')
	streamlines = load_tractogram(streamlines_path,ref_anat)

	# load classification
	print('loading classification')
	classification = sio.loadmat(classification_path)

	# extract names and indices from classification
	names = list(np.ravel(list(classification['classification'][0]['names'][0][0])))
	indices = classification['classification'][0]['index'][0][0]

	# define metrics to use for reorienting streamlines using quickbundles 
	feature = ResampleFeature(nb_points=n_points)
	metric = AveragePointwiseEuclideanMetric(feature)

	# tracts
	tracts = {}
	images_json = {}
	images_json['images'] = []

	# load measures
	df_measures = {}
	for measures in measures_path:
		measure_name = measures.split('/')[-1].split('.nii.gz')[0]
		print('loading measure %s' %measure_name)
		df_measures[measure_name] = nib.load(measures)
		tmp_binary = df_measures[measure_name].get_fdata()[df_measures[measure_name].get_fdata() > 0]
		if np.median(tmp_binary) < 0.01:
			df_measures[measure_name].get_fdata() = df_measures[measure_name].get_fdata() * 1000

	# loop through tracts
	for bundles in range(len(names)):
		print('computing profile for tract %s' %names[bundles])
		
		tracts[names[bundles]] = {}

		# extract streamline data from classificaiton and input streamlines
		tract_index_value = bundles+1
		tract_indices = [ f for f in range(len(indices)) if indices[f] == tract_index_value ]
		fg = streamlines.streamlines[tract_indices]

		# run quickbundles, find centroid, and reorient streamlines
		qb = QuickBundles(np.inf,metric=metric)
		tract_cluster = qb.cluster(fg)
		centroid_cluster = tract_cluster.centroids[0]
		oriented_tract = dts.orient_by_streamline(fg,centroid_cluster)

		# resample to same number of points
		fgarray = dts.set_number_of_points(oriented_tract,n_points)

		# calculate weights
		oriented_tract_weights = dsa.gaussian_weights(fgarray,n_points)

		# loop through measures
		for measures in measures_path:
			measure_name = measures.split('/')[-1].split('.nii.gz')[0]
			print('computing measure %s' %measure_name)

			# get values for each streamline and node
			values = dsa.values_from_volume(df_measures[measure_name].get_fdata(),fgarray,df_measures[measure_name].affine)

			# compute weighted mean
			tracts[names[bundles]][measure_name+'_mean'] = np.sum(oriented_tract_weights * values,0)

			# compute standard deviation
			tracts[names[bundles]][measure_name+'_sd'] = np.std(values,0)

			# output png image of profile for the tract
			print('generating profile images for tract')
			images_json = generateImagesDatatype(images_json,tracts[names[bundles]][measure_name+'_mean'],names[bundles],measure_name,out_path+'/images/%s_%s.png' %(names[bundles],measure_name))
		
		# centroid x,y, and z
		tracts[names[bundles]]['x_coords'] = centroid_cluster[:,0]
		tracts[names[bundles]]['y_coords'] = centroid_cluster[:,1]
		tracts[names[bundles]]['z_coords'] = centroid_cluster[:,2]

		# output profile datatype csv for the tract
		print('generating profiles datatype for tract')
		generateProfilesDatatype(tracts[names[bundles]],out_path+'/profiles/%s_profiles.csv' %names[bundles])		
	
	# dump images.json
	pd.Series(images_json).to_json(out_path+'/images/images.json',orient='index')

	# generate tractmeasures datatype for all tracts
	print('generating tractmeasures datatype')
	generateTractmeasuresDatatype(subjectID,names,list(df_measures.keys()),tracts,out_path+'/tractmeasures/output_FiberStats.csv')

def main():

	import os,sys
	import json

	# load config
	with open('config.json','r') as config_f:
		config = json.load(config_f)

	# make output directories
	if not os.path.exists('./images'):
		os.mkdir('./images')
	if not os.path.exists('./profiles'):
		os.mkdir('./profiles')
	if not os.path.exists('./tractmeasures'):
		os.mkdir('./tractmeasures')

	# define paths and variables
	subjectID = config['_inputs'][0]['meta']['subject']
	streamlines_path = config['track']
	classification_path = config['classification']
	reference_anat_path = config['dwi']
	n_points = config['num_nodes']

	out_path = './'

	# loop through measures to create measures_path
	df_measures = []
	
	# dti
	fa = config['fa']
	if os.path.isfile(fa):
		df_measures = df_measures+['ad','fa','md','rd']

	# dki
	ga = config['ga']
	if os.path.isfile(ga):
		df_measures = df_measures+['ga','ak','mk','rk']

	# noddi
	odi = config['odi']
	if os.path.isfile(odi):
		df_measures = df_measures+['ndi','odi','isovf']

	measure_path = [ config[f] for f in df_measures ]

	# compute tract profiles and output appropriate data
	computeTractProfiles(subjectID,reference_anat_path,streamlines_path,classification_path,measure_path,n_points,out_path)

if __name__ == '__main__':
	main()
