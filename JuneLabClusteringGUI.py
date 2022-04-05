from tkinter import *
from tkinter import ttk
from tkinter import messagebox
from tkinter import filedialog
import tkinter as tk
from tkinter.ttk import Progressbar
import math

from matplotlib.pyplot import text
from scipy.spatial import distance
import GuiBackground as GB
from GUIUtils import GUIUtils as GU 
import logging
import time
import getpass
import os
import fpdf
import webbrowser


class JuneLabClusteringGUI(ttk.Frame):
	def __init__(self, master=None):
		super().__init__(master)
		self.grid(column=0, row=0, sticky=(N, W, E, S))
		self.rowconfigure(0, weight =2)
		self.columnconfigure(1, weight =2)
		self.columnconfigure(2, weight =2)
		self.rowconfigure(1, weight = 2)
		self.rowconfigure(2, weight = 2)
		self.rowconfigure(3, weight = 2)
		self.rowconfigure(4, weight = 2)
		self.columnconfigure(3, weight=2) 
		self.create_widgets()

	def create_widgets(self):
		self.style = ttk.Style()
		self.style.configure("RW.TLabel", foreground="#f03d33",font=("TkHeadingFont",30))
		self.style.configure("RW.TButton", padding=15, borderwidth=15, foreground="black", background="#000000",font=("Arial",14))
		self.JuneLab = ttk.Label(self, text="Welcome to the June Lab Clustering GUI",style="RW.TLabel").grid(column=0,row=0,columnspan=4)
		self.clust = ttk.Button(self,text="Create Clustergram",style="RW.TButton",command=self.createClustergram).grid(column=1,row=1, sticky=(N,S,E,W))
		#Create a button to allow the user to create a medians file for better clustering results. 
		self.med = ttk.Button(self, text="Group Medians", style="RW.TButton", command=self.medians).grid(column=1, row=3, sticky =(N,S,E,W))
		#Create a button to allow the user to compare the four most common linkage functions.
		self.link = ttk.Button(self, text="Compare Linkage Functions", style="RW.TButton", command=self.linkages).grid(column=2, row=1,sticky =(N,S,E,W))
		#Create a button to allow the user to validate the appropriate number of clusters needed for a given set of metabolites.
		self.val = ttk.Button(self, text="Compound Match-Up", style="RW.TButton", command=self.compound).grid(column=3, row=2, sticky =(N,S,E,W))
		#Create a button to allow the user to create the peaks to pathways files needed to analyze the peaks to pathways in Mummichog
		self.peak = ttk.Button(self, text="Peaks to Pathways", style="RW.TButton", command=self.P2P).grid(column=2, row=2, sticky =(N,S,E,W))
		#Create a button to allow the user to check the integrity of the data downloaded from Metaboanalysts Volcano plot results. 
		self.integrity = ttk.Button(self, text="Data Integrity", style="RW.TButton", command=self.integrity).grid(column=2, row=3, sticky =(N,S,E,W))
		#Create a button to allow the user to do an Ensemble clustering on the data.
		self.ensemble = ttk.Button(self,text="Ensemble Clustering", style="RW.TButton", command=self.ensemble).grid(column=1,row=2,sticky=(N,S,E,W))
		#Create a button to allow the user to create a minimum spanning tree on data
		self.mst = ttk.Button(self,text='Minimum Spanning Tree', style="RW.TButton", command=self.mst).grid(column=3,row=1,sticky=(N,S,E,W))
		#Create a button for the generation of a report
		self.generate = ttk.Button(self,text='Generate PDF Report', style="RW.TButton", command=self.generate).grid(column=3,row=3,sticky=(N,S,E,W))
		#Create a button for the users to submit requests. 
		self.request = ttk.Button(self, text="Submit Request", style = "RW.TButton", command=self.userRequest).grid(column=2, row=4, sticky=(N,S,E,W))
		#Create a button for the selection of clusters
		self.selection = ttk.Button(self, text="Cluster Selection",style = "RW.TButton",command=self.clusterSelection).grid(column=1, row=4,sticky=(N,S,E,W))
		#Create a button for locally weighted clustering
		self.localWeight = ttk.Button(self,text="Locally Weighted Ensemble",style="RW.TButton",command=self.localWeighted).grid(column=3, row=4,sticky=(N,S,E,W))
		# pad each widget with 5 pixels on each side to ensure that the buttons do not stay together. 
		for child in self.winfo_children(): child.grid_configure(padx=5, pady=5)


	def home(self):
		#remove old objects and put the home objects back on grid
		objects = self.grid_slaves()
		for i in objects:
			i.grid_remove()
		widgets = self.winfo_children()

		widgetDict = {}
		for i in range(13):
			#create a dictionary of the widgets from home window
			widgetDict[i] = widgets[i]

		widgetDict[0].grid(column=0,row=0,columnspan=4)
		widgetDict[1].grid(column=1,row=1, sticky=(N,S,E,W))
		widgetDict[2].grid(column=1, row=3, sticky =(N,S,E,W))
		widgetDict[3].grid(column=2, row=1,sticky =(N,S,E,W))
		widgetDict[4].grid(column=3, row=2, sticky =(N,S,E,W))
		widgetDict[5].grid(column=2, row=2, sticky =(N,S,E,W))
		widgetDict[6].grid(column=2, row=3, sticky =(N,S,E,W))
		widgetDict[7].grid(column=1,row=2,sticky=(N,S,E,W))
		widgetDict[8].grid(column=3,row=1,sticky=(N,S,E,W))
		widgetDict[9].grid(column=3,row=3,sticky=(N,S,E,W))
		widgetDict[10].grid(column=2, row=4, sticky=(N,S,E,W))
		widgetDict[11].grid(column=1, row=4,sticky=(N,S,E,W))
		widgetDict[12].grid(column=3, row=4, sticky=(N,S,E,W))

		count = -1
		for child in self.winfo_children():
			#add padding to the current widgets
			count += 1
			if count < 13:
				child.grid_configure(padx=5,pady=5)


	def createClustergram(self):
		def linkageOutput(*args):
			#grab the current selection of the list
			global selection
			selection = distListBox.curselection()
			curLink = linkageList[selection[0]]
			if curLink == 'ward':
				lenList = len(sampleListBox.get(0,tk.END))
				if lenList > 0:
					sampleListBox.delete(0,lenList-1)

				sampleListBox.insert(0,distList[0])
				
				self.sampleListBox = sampleListBox
				self.sampleListBox.grid(column=2,row=2,columnspan=1)
			else:
				lenList = len(sampleListBox.get(0,tk.END))
				if lenList > 0:
					sampleListBox.delete(0,lenList-1)

				for i in range(len(distList)):
					sampleListBox.insert(i,distList[i])

				#bind the output back to the GUI.
				self.sampleListBox = sampleListBox
				self.sampleListBox.grid(column=2,row=2,columnspan=1)
			return selection

		def colorMap(*args):
			global selection1
			selection1 = sampleListBox.curselection()
			#create a list of the color map options
			lenList = len(colorListBox.get(0,tk.END))
			if lenList > 0:
				colorListBox.delete(0,lenList-1)

			for i in range(len(colorList)):
				colorListBox.insert(i,colorList[i])

			#bind the output back to the GUI
			self.colorListBox = colorListBox
			self.colorListBox.grid(column=3,row=2,columnspan=1)
			return selection1

		def dataTransform(*args):
			#dataTransform
			global selection2
			selection2 = colorListBox.curselection()

			#put the data transform options into the list
			lenList = len(transformListBox.get(0,tk.END))
			if lenList > 0:
				transformListBox.delete(0,lenList-1)

			for i in range(len(transformList)):
				transformListBox.insert(i,transformList[i])

			self.transformListBox = transformListBox
			self.transformListBox.grid(column=1,row=4,columnspan=1)
			return selection2

		def dataScale(*args):
			#dataScaling
			global selection3
			selection3 = transformListBox.curselection()
			print(selection,selection1,selection2,selection3)
			#put the data scaling optoins into the list
			lenList = len(scaleListBox.get(0,tk.END))
			if lenList > 0:
				scaleListBox.delete(0,lenList-1)

			for i in range(len(scaleList)):
				scaleListBox.insert(i,scaleList[i])

			self.scaleListBox = scaleListBox
			self.scaleListBox.grid(column=2,row=4,columnspan=1)
			return selection3


		def submit(*args):
			#submit the function output to the
			#selection = distListBox.curselection() 
			selection4 = scaleListBox.curselection()
			dist = distList[selection1[0]]
			link = linkageList[selection[0]]
			color = colorList[selection2[0]]
			transform = transformList[selection3[0]]
			scale = scaleList[selection4[0]]
			GU.createClustergram(0,link,dist,color,transform=transform,scale=scale)

		def cmapO(*args):
			#send users to webpage of 
			webbrowser.open('https://matplotlib.org/stable/tutorials/colors/colormaps.html')

		objects = self.grid_slaves()
		for i in objects:
			i.grid_forget()

		#create widgets for the clustergram function input. 
		self.JuneLab = ttk.Label(self, text="Clustergram Input",font=("TkHeadingFont",36)).grid(column=1,row=0,sticky=(N),columnspan=3)
		self.Linkage = ttk.Label(self, text="Linkage",font=("TkHeadingFont",12)).grid(column=1,row=1,sticky=(N))
		self.Distance = ttk.Label(self, text="Distance Measure",font=("TkHeadingFont",12)).grid(column=2,row=1,sticky=(N))
		self.Color = ttk.Label(self,text="Color-Map", font=("TkHeadingFont",12)).grid(column=3,row=1,sticky=(N))
		self.Transform = ttk.Label(self,text="Transform", font=("TkHeadingFont",12)).grid(column=1,row=3,sticky=(N))
		self.Scale = ttk.Label(self,text="Scale",font=('TkHeadingFont',12)).grid(column=2,row=3,sticky=(N))
		self.homepage = ttk.Button(self,text="Return to Home",command=self.home).grid(column=2, row=6,sticky=(N),columnspan=1)
		self.submit = ttk.Button(self,text="Submit", command=submit).grid(column=2, row=5,sticky=(N),columnspan=1)
		self.cmapW = ttk.Button(self,text="ColorMap Options", command=cmapO).grid(column=3,row=3,sticky=(N),columnspan=1)
		distListBox = Listbox(self,height=8)
		sampleListBox = Listbox(self,height=8)
		colorListBox = Listbox(self,height=8)
		transformListBox = Listbox(self, height=8)
		scaleListBox = Listbox(self, height=8)
		
		#Create the lists of available options for selection 
		linkageList = ('single','ward','complete','average')
		distList = ('euclidean','seuclidean','sqeuclidean','cosine','chebyshev','correlation','canberra','braycurtis','minkowski','cityblock')
		colorList = ('viridis', 'plasma', 'inferno', 'magma', 'cividis','Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                      'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                      'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn','Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                      'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                      'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn','PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu',
                      'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic','twilight', 'twilight_shifted', 'hsv')
		transformList = ('None','Log transformation', 'Square root transformation', 'Cube root transformation')
		scaleList = ('None', 'Mean centering', 'Auto Scaling', 'Pareto Scaling', 'Range Scaling')

		
		linkNames = StringVar(value=linkageList)
		distNames = StringVar(value=distList)
		colorNames = StringVar(value=colorList)
		transformNames = StringVar(value=transformList)
		scaleNames = StringVar(value=scaleList)

		
		#input the linkage function values into the box
		for i in range(len(linkageList)):
			distListBox.insert(i,linkageList[i])

		distListBox.bind('<Double-1>',linkageOutput)
		sampleListBox.bind('<Double-1>',colorMap)
		colorListBox.bind('<Double-1>',dataTransform)
		transformListBox.bind('<Double-1>', dataScale)
		self.distListBox = distListBox
		self.distListBox.grid(column=1,row=2,columnspan=1)
		self.sampleListBox = sampleListBox
		self.sampleListBox.grid(column=2,row=2,columnspan=1)
		self.colorListBox = colorListBox
		self.colorListBox.grid(column=3,row=2,columnspan=1)
		self.transformListBox = transformListBox
		self.transformListBox.grid(column=1,row=4,columnspan=1)
		self.scaleListBox = scaleListBox
		self.scaleListBox.grid(column=2,row=4,columnspan=1)

	def clusterSelection(self):
		def linkageOutput(*args):
			#grab the current selection of the list
			global selection
			selection = distListBox.curselection()
			curLink = linkageList[selection[0]]
			if curLink == 'ward':
				lenList = len(sampleListBox.get(0,tk.END))
				if lenList > 0:
					sampleListBox.delete(0,lenList-1)

				sampleListBox.insert(0,distList[0])
				
				self.sampleListBox = sampleListBox
				self.sampleListBox.grid(column=2,row=2,columnspan=1)
			else:
				lenList = len(sampleListBox.get(0,tk.END))
				if lenList > 0:
					sampleListBox.delete(0,lenList-1)

				for i in range(len(distList)):
					sampleListBox.insert(i,distList[i])

				#bind the output back to the GUI.
				self.sampleListBox = sampleListBox
				self.sampleListBox.grid(column=2,row=2,columnspan=1)
			return selection

		def colorMap(*args):
			global selection1
			selection1 = sampleListBox.curselection()
			#create a list of the color map options
			lenList = len(colorListBox.get(0,tk.END))
			if lenList > 0:
				colorListBox.delete(0,lenList-1)

			for i in range(len(colorList)):
				colorListBox.insert(i,colorList[i])

			#bind the output back to the GUI
			self.colorListBox = colorListBox
			self.colorListBox.grid(column=3,row=2,columnspan=1)
			return selection1

		def dataTransform(*args):
			#dataTransform
			global selection2
			selection2 = colorListBox.curselection()

			#put the data transform options into the list
			lenList = len(transformListBox.get(0,tk.END))
			if lenList > 0:
				transformListBox.delete(0,lenList-1)

			for i in range(len(transformList)):
				transformListBox.insert(i,transformList[i])

			self.transformListBox = transformListBox
			self.transformListBox.grid(column=1,row=4,columnspan=1)
			return selection2

		def dataScale(*args):
			#dataScaling
			global selection3
			selection3 = transformListBox.curselection()
			print(selection,selection1,selection2,selection3)
			#put the data scaling optoins into the list
			lenList = len(scaleListBox.get(0,tk.END))
			if lenList > 0:
				scaleListBox.delete(0,lenList-1)

			for i in range(len(scaleList)):
				scaleListBox.insert(i,scaleList[i])

			self.scaleListBox = scaleListBox
			self.scaleListBox.grid(column=2,row=4,columnspan=1)
			return selection3


		def submit(*args):
			#submit the function output to the
			#selection = distListBox.curselection() 
			selection4 = scaleListBox.curselection()
			dist = distList[selection1[0]]
			link = linkageList[selection[0]]
			color = colorList[selection2[0]]
			transform = transformList[selection3[0]]
			scale = scaleList[selection4[0]]
			GU.selectClusters(link,dist,transform=transform, scale=scale,cmap=color)

		def cmapO(*args):
			#send users to webpage of 
			webbrowser.open('https://matplotlib.org/stable/tutorials/colors/colormaps.html')

		objects = self.grid_slaves()
		for i in objects:
			i.grid_forget()

		#create widgets for the clustergram function input. 
		self.JuneLab = ttk.Label(self, text="Clustergram Input",font=("TkHeadingFont",36)).grid(column=1,row=0,sticky=(N),columnspan=3)
		self.Linkage = ttk.Label(self, text="Linkage",font=("TkHeadingFont",12)).grid(column=1,row=1,sticky=(N))
		self.Distance = ttk.Label(self, text="Distance Measure",font=("TkHeadingFont",12)).grid(column=2,row=1,sticky=(N))
		self.Color = ttk.Label(self,text="Color-Map", font=("TkHeadingFont",12)).grid(column=3,row=1,sticky=(N))
		self.Transform = ttk.Label(self,text="Transform", font=("TkHeadingFont",12)).grid(column=1,row=3,sticky=(N))
		self.Scale = ttk.Label(self,text="Scale",font=('TkHeadingFont',12)).grid(column=2,row=3,sticky=(N))
		self.homepage = ttk.Button(self,text="Return to Home",command=self.home).grid(column=2, row=6,sticky=(N),columnspan=1)
		self.submit = ttk.Button(self,text="Submit", command=submit).grid(column=2, row=5,sticky=(N),columnspan=1)
		self.cmapW = ttk.Button(self,text="ColorMap Options", command=cmapO).grid(column=3,row=3,sticky=(N),columnspan=1)
		distListBox = Listbox(self,height=8)
		sampleListBox = Listbox(self,height=8)
		colorListBox = Listbox(self,height=8)
		transformListBox = Listbox(self, height=8)
		scaleListBox = Listbox(self, height=8)
		
		#Create the lists of available options for selection 
		linkageList = ('single','ward','complete','average')
		distList = ('euclidean','seuclidean','sqeuclidean','cosine','chebyshev','correlation','canberra','braycurtis','minkowski','cityblock')
		colorList = ('viridis', 'plasma', 'inferno', 'magma', 'cividis','Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                      'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                      'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn','Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                      'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                      'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn','PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu',
                      'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic','twilight', 'twilight_shifted', 'hsv')
		transformList = ('None','Log transformation', 'Square root transformation', 'Cube root transformation')
		scaleList = ('None', 'Mean centering', 'Auto Scaling', 'Pareto Scaling', 'Range Scaling')

		
		linkNames = StringVar(value=linkageList)
		distNames = StringVar(value=distList)
		colorNames = StringVar(value=colorList)
		transformNames = StringVar(value=transformList)
		scaleNames = StringVar(value=scaleList)

		
		#input the linkage function values into the box
		for i in range(len(linkageList)):
			distListBox.insert(i,linkageList[i])

		distListBox.bind('<Double-1>',linkageOutput)
		sampleListBox.bind('<Double-1>',colorMap)
		colorListBox.bind('<Double-1>',dataTransform)
		transformListBox.bind('<Double-1>', dataScale)
		self.distListBox = distListBox
		self.distListBox.grid(column=1,row=2,columnspan=1)
		self.sampleListBox = sampleListBox
		self.sampleListBox.grid(column=2,row=2,columnspan=1)
		self.colorListBox = colorListBox
		self.colorListBox.grid(column=3,row=2,columnspan=1)
		self.transformListBox = transformListBox
		self.transformListBox.grid(column=1,row=4,columnspan=1)
		self.scaleListBox = scaleListBox
		self.scaleListBox.grid(column=2,row=4,columnspan=1)

	def medians(self):
		global rmZeros
		#function to 
		def groupMedians(*args):
			rmZeros = var1.get()
			GU.groupMedians(rmZeros=rmZeros)
		#get rid of the objects in GUI window.
		objects = self.grid_slaves()
		for i in objects:
			i.grid_forget()

		#create widgets for the group medians.
		self.GroupLab = ttk.Label(self, text="Group Medians", font=("TkHeadingFont",36)).grid(column=1,row=0,sticky=(N),columnspan=2)
		self.homeGroup = ttk.Button(self,text="Return to Home",command=self.home).grid(column=1,row=4,sticky=(N),columnspan=2)
		var1 = IntVar()
		self.zeroRemove = ttk.Checkbutton(self,text="Remove Zeros?",variable=var1).grid(column=1,row=2,sticky=(N),columnspan=2)

		self.GM = ttk.Button(self,text="Select file",command=groupMedians).grid(column=1,row=3,sticky=(N),columnspan=2)


	def linkages(self):
		def distFunc(*args):
			#make a global distance variable
			global distance
			distance = self.dist.get()
			
			#given the selection of a distance measure how many comparisons are possible. 
			if distance == 'euclidean':
				#give the full list to the second combobox
				values = [1,2,3,4]
				num_comps = StringVar()
				self.numComps = ttk.Combobox(self,values=values,textvariable=num_comps)
			else:
				#give only three values to the second combobox
				values = [1,2,3]
				num_comps = StringVar()
				self.numComps = ttk.Combobox(self,values=values,textvariable=num_comps)
			
			#place the combobox on the GUI
			self.numComps.bind('<<ComboboxSelected>>', numCompsFunc)
			self.numComps.grid(column=2,row=2)
		
		def numCompsFunc(*args):
			#make a numberComps global variable
			global numberComps
			numberComps = self.numComps.get()

			#create a list of linkage options
			linkOpts = ['ward-single','ward-complete','ward-average','single-complete','single-average','complete-average',\
						'ward-single-complete','ward-single-average','ward-complete-average','single-complete-average',\
						'ward-single-complete-average']

			if numberComps == '1':
				#give the linkage funtions combobox a list of the linkage functions. 
				linkages = ['ward','single','complete','average']

				#check for the euclidean measure again.
				if distance == 'euclidean':
					#give the values of linkages[1:3]
					value = linkages
					linkage = StringVar()
					self.linkage = ttk.Combobox(self,values=value,textvariable=linkage)
				else:
					value = linkages[3]
					linkage = StringVar()
					self.linkage = ttk.Combobox(self,values = value,textvariable=linkage)

				#bind and place the combobox on the GUI
				self.linkage.bind('<<ComboboxSelected>>', linkageComp)
				self.linkage.grid(column=3,row=2)

			elif numberComps == '2':
				#give the linkage functions combobox a subset based upon distance measure. 
				if distance == 'euclidean':
					#give the list values from 0 to 5
					value = linkOpts[0:6]
					linkage = StringVar()
					self.linkage = ttk.Combobox(self,values=value,textvariable=linkage)

				else:
					#give the list values from 3 to 5
					value = linkOpts[3:5]
					linkage = StringVar()
					self.linkage = ttk.Combobox(self,values=value,textvariable=linkage)
				
				#bind and place the combobox on the GUI
				self.linkage.bind('<<ComboboxSelected>>', linkageComp)
				self.linkage.grid(column=3,row=2)

			elif numberComps == '3':
				#give the linkage functions combobox a subset based upon distance measure. 
				if distance == 'euclidean':
					#give the list values from 0 to 5
					value = linkOpts[6:10]
					linkage = StringVar()
					self.linkage = ttk.Combobox(self,values=value,textvariable=linkage)

				else:
					#give the list values from 3 to 5
					value = linkOpts[9]
					linkage = StringVar()
					self.linkage = ttk.Combobox(self,values=value,textvariable=linkage)
				
				#bind and place the combobox on the GUI
				self.linkage.bind('<<ComboboxSelected>>', linkageComp)
				self.linkage.grid(column=3,row=2)

			elif numberComps == '4':
				#give the list values from 3 to 5
				value = linkOpts[10]
				linkage = StringVar()
				self.linkage = ttk.Combobox(self,values=value,textvariable=linkage)
				
				#bind and place the combobox on the GUI
				self.linkage.bind('<<ComboboxSelected>>', linkageComp)
				self.linkage.grid(column=3,row=2)

		def linkageComp(*args):
			#base linkage list for 4 comparisons.
			linkageList = ['ward','single','complete','average']
			
			#create an empty list for linkage functions
			linkList = []

			#count the number of dashes in the string from combobox
			selection = self.linkage.get()
			locs = []
			for i in range(len(selection)):
				#find the dash locations
				if selection[i] == '-':
					locs.append(i)

			if len(locs) > 0:
				firstLetter = 0
				for i in range(len(selection)):
					#check vthe current string value for '-'
					if selection[i] == '-':
						lastLetter = i
						#get current linkage
						curLink = selection[firstLetter:lastLetter]
						linkList.append(curLink)
						firstLetter = i+1
					elif i == len(selection)-1:
						lastLetter = i+1
						curLink = selection[firstLetter:lastLetter]
						linkList.append(curLink)

			else:
				#append the selection to the linkage list
				linkList.append(selection)
			
			# #send the parameters for linkage comparison 	
			file = filedialog.askopenfilename()	
			GU.linkageComparison(file,numberComps,linkList,distance)



		objects = self.grid_slaves()
		for i in objects:
			i.grid_forget()

		#create a list of values from 1 to 4
		numLinkComps = [1,2,3,4]

		#create widgets for the clustergram function input. 
		self.JuneLab = ttk.Label(self, text="Linkage Comparison",font=("TkHeadingFont",36)).grid(column=1,row=0,sticky=(N),columnspan=3)
		self.numCompsLab = ttk.Label(self, text="Number of comparisons",font=("TkHeadingFont",12)).grid(column=2,row=1)
		self.distLab = ttk.Label(self, text="Distance measure",font=("TkHeadingFont",12)).grid(column=1,row=1)
		self.linkLab = ttk.Label(self,text="Linkage functions",font=("TkHeadingFont",12)).grid(column=3,row=1)
		self.home1 = ttk.Button(self,text="Return to Home",command=self.home).grid(column=1,row=4, sticky=(N),columnspan=3)

		linkages = StringVar()
		#create the distance measure combobox first, then update the GUI as the user selects the distance, measure than number of linkage comps. 
		distances = StringVar()
		distList = ('euclidean','seuclidean','sqeuclidean','cosine','chebyshev','correlation','canberra','braycurtis','minkowski','cityblock')
		self.dist = ttk.Combobox(self,values = distList,textvariable=distances)
		
		#bind the combobox for distance measures to the selection of distance measure. 
		self.dist.bind('<<ComboboxSelected>>', distFunc)
		self.dist.grid(column=1,row=2)

	def compound(self):
	    #ask the user to select a clustergram file to run through a validition study.
	    #Waiting on confirmation...
		def allCompounds(*args):
			GU.compoundMatchUp(typeFile='all')

		def enrichment(*args):
			GU.compoundMatchUp(typeFile='enrich')

		objects = self.grid_slaves()
		for i in objects:
			i.grid_forget()

		self.style.configure("RW.TButton", padding=15, borderwidth=15, foreground="black", background="#000000",font=("Arial",14))
		self.JuneLab = ttk.Label(self, text ="Compound Match-Up",font=("TkHeadingFont",36)).grid(column=1,row=0,sticky=(N),columnspan=3)
		self.allCompounds = ttk.Button(self,text="All Compounds",command=allCompounds,style="RW.TButton").grid(column=1,row=1, sticky=(N),columnspan=3)
		self.enrichmentCompounds = ttk.Button(self,text="Enrichment Compounds",command = enrichment,style="RW.TButton").grid(column=1,row=2, sticky=(N),columnspan=3)
		self.home1 = ttk.Button(self,text="Return to Home",command=self.home,style="RW.TButton").grid(column=1,row=3, sticky=(N),columnspan=3)


	def P2P(self):
	    #Waiting until the MST functionality is complete. 
	    GU.peaksToPathways()

	def integrity(self):
	    #ask the user to select a volcano plot file to check the integrity of the data against. 
	    filename = filedialog.askopenfilename()
	    GU.dataIntegrity(filename)

	def mstF(self):
		numClust = GU.MST(func='ensemble')
		numClust = int(numClust)
		print(numClust)
		GU.ensembleClustering(optNum = numClust)

	def ensemble(self):
		def ensembleOutput(*args):
			#grab the current selection of the list
			global selection
			selection = optClustBox.curselection()
			selection = int(selection[0])+1

		def ensembleMin(*args):
			#grab the selection and minimum number of clusters
			global selection1
			selection1 = minNumMetabClust.curselection()
			selection1 = int(selection1[0])

		def ensembleGo(*args):
			selection2 = colorMapEnsem.curselection()
			selection2 = colorList[selection2[0]]
			print(selection)
			print(selection1)
			print(selection2)
			#send to ensemble clustering algorithm
			GU.ensembleClustering(optNum= selection, minMetabs= selection1, colorMap = selection2)

		objects = self.grid_slaves()
		for i in objects:
			i.grid_forget()

		colorList = ('viridis', 'plasma', 'inferno', 'magma', 'cividis','Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                      'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                      'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn','Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                      'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                      'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn','PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu',
                      'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic','twilight', 'twilight_shifted', 'hsv')	
		#create widgets for the clustergram function input. 
		self.JuneLab = ttk.Label(self, text="Number of Clusters",font=("TkHeadingFont",18)).grid(column=0,row=0,sticky=(N))
		self.MetabNumPClustL = ttk.Label(self, text="Min. number of metabolites (per cluster)?", font=("TkHeadingFont",18)).grid(column=0,row=3,sticky=(N))
		self.colMap = ttk.Label(self, text="ColorMap", font=("TkHeadingFont",18)).grid(column=0,row=5,sticky=(N))
		self.mstF1 = ttk.Button(self,text="Run MST",command=self.mstF).grid(column=0,row=7,sticky=(N),columnspan=2)
		self.home1 = ttk.Button(self,text="Return to Home",command=self.home).grid(column=0,row=8, sticky=(N),columnspan=2)

		#create list box of the ensemble optimal clusters
		optClustBox = Listbox(self,height=5)
		minNumMetabClust = Listbox(self, height=5)
		colorMapEnsem = Listbox(self, height=5)

		#Create the lists of available options for selection 
		optClusters = tuple(range(1,101))
		minMetab = tuple(range(0,101))
		distNames = StringVar(value=optClusters)
		minNum = StringVar(value=minMetab)
		colMap = StringVar(value=colorList)


		for i in range(len(optClusters)):
			optClustBox.insert(i,optClusters[i])

		for i in range(len(minMetab)):
			minNumMetabClust.insert(i,minMetab[i])

		for i in range(len(colorList)):
			colorMapEnsem.insert(i,colorList[i])

		#create binding event and 
		optClustBox.bind('<Double-1>',ensembleOutput)
		minNumMetabClust.bind('<Double-1>',ensembleMin)
		colorMapEnsem.bind('<Double-1>',ensembleGo)
		self.optClustBox = optClustBox
		self.optClustBox.grid(column=0,row=1,columnspan=1)
		self.minNumMetabClust = minNumMetabClust
		self.minNumMetabClust.grid(column=0,row=4,columnspan=1)
		self.colorMapEnsem = colorMapEnsem
		self.colorMapEnsem.grid(column=0,row=6,columnspan=1)


	def mst(self):
		def valType(*args):
			selection = valTypeBox.curselection()
			index = selection[0]
			#print(type(index))
			selection =valList[selection[0]]
			print(selection)
			#send the user to the minimum spanning tree function.
			GU.MST(self,func='base')
		

		#delete previous objects on the GUI
		objects = self.grid_slaves()
		for i in objects:
			i.grid_forget()

		#create widgets for the clustergram function input. 
		self.JuneLab = ttk.Label(self, text="Validation Measure",font=("TkHeadingFont",18)).grid(column=0,row=0,sticky=(N))
		self.home1 = ttk.Button(self,text="Return to Home",command=self.home).grid(column=0,row=3, sticky=(N),columnspan=2)

		#validation index list (MST-based, DBI)
		valList = ('MST-based','DBI')
		valTypeBox = Listbox(self,height=5)

		valName = StringVar(value=valTypeBox)

		for i in range(len(valList)):
			valTypeBox.insert(i, valList[i])

		#create a binding event
		valTypeBox.bind('<Double-1>',valType)
		self.valTypeBox = valTypeBox
		self.valTypeBox.grid(column=0,row=1,sticky=(N))

		
		


	def generate(self):
	    #ask the user to select the file in which will be used to create a minimum spanning tree.
	    GU.PDFGenerator()

	def localWeighted(self):
		#send the user to locally weighted ensemble clustering
		GU.localWeighted()
	
	def userRequest(self):
		def generateRequest(*args):
			#create the pdf and title for each page
			pdf = fpdf.FPDF('P','mm','Letter')
			title = self.title.get()
			body = self.body.get()

			#Create the title and set the default font
			directory = 'C:/Users/Public/Documents/Requests'
			os.chdir(directory)

			#determine the current user
			curUser = getpass.getuser()
			curUser = GB.who(curUser)
			
			#create first page
			pdf.add_page()
			title += '-' + curUser
			pdf.set_font('Arial', 'B', 24)
			pdf.cell(197, 10, title, 0, 0, 'C')
			pdf.line(5,20,200,20)
			pdf.ln(15)
			pdf.set_font('Arial',style='',size=18)

			pdf.multi_cell(197,8,body,0,0,'J')
			ending = '.pdf'
			fileName = ''
			curTime = time.strftime("%a_%b_%d_%Y_%H_%M")
			fileName += curUser + '_' + curTime + ending
			pdf.output(fileName, 'F')
			logging.info(': Sucessfully created a pdf of the results!')
			logging.info(': Leaving the pdf User Request Function!')

		objects = self.grid_slaves()
		for i in objects:
			i.grid_forget()

		self.title = tk.StringVar()
		self.body = tk.StringVar()
		self.entrytitle = ttk.Entry(self,textvariable=self.title).grid(column=0,row=1,columnspan=4, pady=3)
		self.entrybody = ttk.Entry(self,textvariable=self.body).grid(column=0,row=3,columnspan=4, pady=3)
		self.labelTitle = ttk.Label(self, text='Title').grid(row=0)
		self.labelBody = ttk.Label(self, text='Description').grid(row=2)
		self.home4 = ttk.Button(self,text="Return to Home",command=self.home).grid(column=0,row=4, sticky=(N),columnspan=2)

		self.buttonSubmit = ttk.Button(self,text = 'Submit Request', command = generateRequest).grid(column=0, columnspan=4)
		

curUser = getpass.getuser()
curUser = GB.who(curUser)
log_time = time.strftime("%a_%b_%d_%Y_%H_%M_%S")
directory = "C:/Users/" + getpass.getuser() + '/Desktop/LogPOutputFiles'
os.chdir(directory)
log_file = curUser + '_' + str(log_time) + '.log' 
logging.basicConfig(filename=log_file,format='%(asctime)s %(message)s',datefmt='%m/%d/%Y %I:%M:%S %p',level=logging.INFO)
logging.info('Started GUI')

if __name__ == '__main__':
	root = tk.Tk()
	app = JuneLabClusteringGUI(master=None)
	app.master.minsize(200,200)
	height = app.master.winfo_screenheight() - 400
	width = (app.master.winfo_screenwidth())*0.9
	width = int(width)
	app.master.maxsize(800,800)
	width = str(width)
	height = str(height)
	screenSize = width +'x' + height + '+0+100'
	app.master.geometry("")
	app.master.resizable(0,0)

	app.master.title("Home")
	app.mainloop()