import logging, time, os, multiprocessing, webbrowser, matplotlib, config
from collections import defaultdict
matplotlib.use("TkAgg")
import tkinter as tk
from tkinter import *
from tkinter import ttk, messagebox, filedialog
from GUIUtils import GUIUtils as GU
import GuiBackground as GB

class CreateToolTip:
    '''This functionality provides a pop-up window to each button explaining the functionality'''
    def __init__(self, widget, text=''):
        self.widget = widget
        self.text   = text
        self.id     = None
        self.tw     = None
        widget.bind("<Enter>",       self.enter)
        widget.bind("<Leave>",       self.leave)
        widget.bind("<ButtonPress>", self.leave)

    def enter(self, e=None):  self._schedule()
    def leave(self, e=None):  self._unschedule(); self._hide()

    def _schedule(self):
        self._unschedule()
        self.id = self.widget.after(500, self._show)

    def _unschedule(self):
        id_, self.id = self.id, None
        if id_: self.widget.after_cancel(id_)

    def _show(self, e=None):
        # x, y, *_ = self.widget.bbox("insert")
        # x += self.widget.winfo_rootx() + 25
        # y += self.widget.winfo_rooty() + 20
        x = self.widget.winfo_rootx() + self.widget.winfo_width() // 2
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 4
        self.tw = tk.Toplevel(self.widget)
        self.tw.wm_overrideredirect(True)
        self.tw.wm_geometry(f"+{x}+{y}")
        tk.Label(self.tw, text=self.text, justify='left',
                 background="#ffffff", relief='solid', borderwidth=1,
                 wraplength=200).pack(ipadx=1)

    def _hide(self):
        tw, self.tw = self.tw, None
        if tw: tw.destroy()

class JuneLabClusteringGUI(ttk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        # Fill the whole window
        self.grid(column=0, row=0, sticky=(N, W, E, S))
        master.columnconfigure(0, weight=1)
        master.rowconfigure(0, weight=1)
        self.startUpPage()

    #  helpers 
    def _clear(self):
        """Destroy every child widget so the frame is completely empty."""
        for w in self.winfo_children():
            w.destroy()

    def _run_with_busy(self, status_text, work_fn, status_row=8):
        """Show a visible busy state while running a long callback."""
        status_lbl = ttk.Label(self, text=status_text, font=("TkHeadingFont", 11))
        status_lbl.grid(column=1, row=status_row, columnspan=3, pady=4)
        self.configure(cursor="watch")
        if self.master is not None:
            self.master.configure(cursor="watch")
        self.update_idletasks()
        try:
            return work_fn()
        finally:
            status_lbl.destroy()
            self.configure(cursor="")
            if self.master is not None:
                self.master.configure(cursor="")
            self.update_idletasks()

    def _configure_home_grid(self):
        """
        Set column/row weights for the main button grid.
        Called once after _clear() before re-populating with buttons.
        """
        # reset any previous config
        for c in range(5):
            self.columnconfigure(c, weight=0)
        for r in range(12):
            self.rowconfigure(r, weight=0)

        # Button grid uses columns 1..4 (4 equal-width button columns).
        # Column 0 acts as the only left padding column.
        self.columnconfigure(0, weight=1)   # left pad
        self.columnconfigure(1, weight=3)
        self.columnconfigure(2, weight=3)
        self.columnconfigure(3, weight=3)
        self.columnconfigure(4, weight=3)   # right button column

        self.rowconfigure(0, weight=1)      # title row
        for r in range(1, 9):
            self.rowconfigure(r, weight=2)  # button rows

    # Sub-page layout: content uses rows 1..(FOOTER_ROW-2); spacer grows;
    # Return to Home sits on FOOTER_ROW (see _sub_page).
    _SUB_SPACER_ROW = 15
    _SUB_FOOTER_ROW = 16

    def _configure_sub_grid(self):
        """Weights for sub-pages (centred layout; footer row set in _sub_page)."""
        for c in range(5):
            self.columnconfigure(c, weight=0)

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=3)
        self.columnconfigure(2, weight=3)
        self.columnconfigure(3, weight=3)
        self.columnconfigure(4, weight=1)

        for r in range(24):
            self.rowconfigure(r, weight=0)

    def _configure_startup_grid(self):
        """Setup screen only: same columns as sub-pages, rows get weight for spacing."""
        for c in range(5):
            self.columnconfigure(c, weight=0)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=3)
        self.columnconfigure(2, weight=3)
        self.columnconfigure(3, weight=3)
        self.columnconfigure(4, weight=1)
        for r in range(24):
            self.rowconfigure(r, weight=0)
        for r in range(1, 9):
            self.rowconfigure(r, weight=1)

    # Wiki links (match Ensemble Clustering labels/URLs)
    _URL_COLORMAP = 'https://matplotlib.org/stable/tutorials/colors/colormaps.html'
    _URL_TRANSFORM = (
        'https://github.com/hisl6802/Transformation-and-Scaling/wiki/Transformations'
    )
    _URL_SCALE = 'https://github.com/hisl6802/Transformation-and-Scaling/wiki/Scaling'

    def _wiki_buttons_color_transform_scale(self, row, start_col=1):
        for col, (label, url) in enumerate(
            [
                ('ColorMap Options', self._URL_COLORMAP),
                ('Transform Options', self._URL_TRANSFORM),
                ('Scale Options', self._URL_SCALE),
            ],
            start=start_col,
        ):
            ttk.Button(
                self, text=label,
                command=lambda u=url: webbrowser.open(u),
            ).grid(column=col, row=row, padx=5, pady=5, sticky=(E, W))

    def _wiki_buttons_transform_scale(self, row, start_col=1):
        for col, (label, url) in enumerate(
            [
                ('Transform Options', self._URL_TRANSFORM),
                ('Scale Options', self._URL_SCALE),
            ],
            start=start_col,
        ):
            ttk.Button(
                self, text=label,
                command=lambda u=url: webbrowser.open(u),
            ).grid(column=col, row=row, padx=5, pady=5, sticky=(E, W))

    def _wiki_buttons_cmap_trans_scale_under_menus(self, row):
        """Wiki under Color-Map / Transform / Scale listbox columns 3–5."""
        for col, (label, url) in (
            (3, ('ColorMap Options', self._URL_COLORMAP)),
            (4, ('Transform Options', self._URL_TRANSFORM)),
            (5, ('Scale Options', self._URL_SCALE)),
        ):
            ttk.Button(
                self, text=label,
                command=lambda u=url: webbrowser.open(u),
            ).grid(column=col, row=row, padx=5, pady=5, sticky=(E, W))

    #  start-up page
    def startUpPage(self):
        self._clear()
        self._configure_startup_grid()

        self.style = ttk.Style()
        self.style.configure("RW.TLabel",  foreground="#000000",
                             font=("TkHeadingFont", 30))
        self.style.configure("RW.TButton", padding=15, borderwidth=15,
                             foreground="gray", background="#000000",
                             font=("Arial", 14))

        numThreads = int(multiprocessing.cpu_count())

        ttk.Label(self, text="Clustering Toolbox Set-Up", style="RW.TLabel").grid(
            column=1, row=0, columnspan=3, pady=20)
        ttk.Label(self, text="Please input your name or a Project name:",
                  font=("TkHeadingFont", 16)).grid(
            column=1, row=1, columnspan=3, sticky=(E, W), pady=2)

        self.name    = tk.StringVar()
        self.threads = tk.StringVar()
        self.output_dir = tk.StringVar()

        ttk.Entry(self, textvariable=self.name).grid(
            column=1, row=2, columnspan=3, sticky=(E, W), padx=20)
        ttk.Label(self, text="Number of threads:",
                  font=("TkHeadingFont", 16)).grid(
            column=1, row=3, columnspan=3, sticky=(E, W), pady=2)
        ttk.Entry(self, textvariable=self.threads).grid(
            column=1, row=4, columnspan=3, sticky=(E, W), padx=20)
        ttk.Label(self,
                  text=f"You have {numThreads} available"
                       "(Using half or less is recommended)",
                  font=("TkHeadingFont", 14)).grid(
            column=1, row=5, columnspan=3, sticky=(E, W), pady=2)

        ttk.Label(self, text="Output folder (logs + results):",
                  font=("TkHeadingFont", 16)).grid(
            column=1, row=6, columnspan=3, sticky=(E, W), pady=(10, 2))
        ttk.Entry(self, textvariable=self.output_dir).grid(
            column=1, row=7, columnspan=2, sticky=(E, W), padx=(20, 5))
        ttk.Button(self, text="Browse…",
                   command=lambda: self.output_dir.set(
                       filedialog.askdirectory() or self.output_dir.get()
                   )).grid(column=3, row=7, sticky=(E, W), padx=(5, 20))

        ttk.Button(self, text="Get Started!",
                   command=self.create_widgets).grid(
            column=1, row=8, columnspan=3, sticky=(E, W), padx=30, pady=(12, 2))

    #  main button page
    def create_widgets(self):
        '''
        This creates the home page full of all of the buttons for the Clustering Toolbox.
        '''
        name       = self.name.get()
        numThreads = self.threads.get()
        out_dir    = (self.output_dir.get() or "").strip()
        try:
            if int(numThreads) <= multiprocessing.cpu_count():
                config.numThreads = int(numThreads)
        except Exception:
            config.numThreads = 2
            messagebox.showinfo(
                title="Number of Threads",
                message="You have been assigned 2 threads.")
        config.name = name
        log_time = time.strftime("%a_%b_%d_%Y_%H_%M_%S")
        if out_dir:
            try:
                os.chdir(out_dir)
                session_dir = os.path.join(os.getcwd(), f"{name}_{log_time}")
                os.makedirs(session_dir, exist_ok=True)
                os.chdir(session_dir)
            except Exception:
                messagebox.showwarning(
                    title="Output folder",
                    message=("Unable to switch to the selected output folder.\n\n"
                             f"Folder: {out_dir}\n\n"
                             "Logs/results will be saved to the current working directory.")
                )

        log_file  = config.name + '_' + log_time + '.log'
        logging.basicConfig(filename=log_file,
                            format='%(asctime)s %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            level=logging.INFO)
        logging.info('Clustering Toolbox appropriately started...')

        self._show_home()

    def _show_home(self):
        """
        Destroy everything, re-apply home grid weights, then rebuild
        all buttons from scratch.  Called by create_widgets AND home().
        """
        self._clear()
        self._configure_home_grid()

        self.style = ttk.Style()
        self.style.configure('RW.TLabel',
                             foreground='#000000',
                             font=('Arial', 28, 'bold'),
                             padding=10)
        self.style.configure('RW.TButton',
                             font=('Arial', 11, 'bold'),
                             foreground='#000000',
                             background='#95a8ba',
                             borderwidth=5,
                             padding=8,
                             relief='raised')

        ttk.Label(self, text="Clustering Toolbox",
                  style="RW.TLabel").grid(
            column=0, row=0, columnspan=5, pady=10)

        # button definitions: (label, command, col, row, tooltip)
        buttons = [
            ("Remove Duplicate Features", self.preprocess,          1, 1,
             "Remove rows with duplicate values in the first column"),
            
            ("Data Integrity",            self.integrity,           2, 1,
             "Handles multiple decimals in first column common in MetaboAnalyst Volcano plot results downloads"),
            
            ("Check Normality",           self.normalityC,          3, 1,
             "Check normality - shows before and after data transformation"),

            ("Compare Linkage Functions", self.linkages,            4, 1,
             "Compare linkage functions used for clustering outputs"),

            ("Create Clustergram",        self.createClustergram,  1, 2,
             "Generate a clustergram heatmap from data"),

            ("Validation Clusters",       self.mst,                2, 2,
             "Run cluster-count validation across 5 metrics (k-means based, DBI, Dunn, PBM, Silhouette, or All). Saves CSV and plots to help pick K."),

            ("Monoclustering Validation", self.monoVal,             3, 2,
             "Validate mono-clustering solutions"),

            ("Group Medians",             self.medians,             4, 2,
             "Calculate group medians from raw data"),

            ("Ensemble Clustering",       self.ensemble,            1, 3,
             "Run ECCO"),

            ("Peaks to Pathways",         self.P2P,                 2, 3,
             "Creates a Peaks to Pathways compatible file for each cluster"),

            ("Compound ID",               self.compound,            3, 3,
             "Lookup compounds from Mummichog output. Select folder of mummichog_matched_compound_*.csv (run MummiBot on P2P files first)."),

            ("CoOcc Heatmap",             self.cooccHeatmap,       4, 3,
             "Create a heatmap from the CoOcc output on the group medians output. Better for smaller datasets."),

            ("Selected Clusters Figure",  self.genSelClustFig,      1, 4,
             "Makes a figure for specific clusters"),

            ("Cluster Selection",         self.clusterSelection,    2, 4,
             "Select clusters interactively"),

            ("Heatmap Analyses",          self.heatmapAnalyses,     3, 4,
             "Make a relative-intensity heatmap using ECCO-identified clusters"),

            ("Build ANOVA Heatmap",       self.anovaHeatMap,        4, 4,
             "Make a heatmap with ECCO results"),

            ("CIs for t-tests",           self.CIsTtest,            1, 5,
             "Create confidence intervals from t-test data"),

            ("MST Optimization",         self.mstF,                 2, 5,
             "Run MST validation (k-means based), suggest optimal K, then run ECCO ensemble with that K"),

            ("Enzyme Look Up",            self.enzymeLookUp,        3, 5,
             "Look up enzymes for found compounds"),

            ("Metaboanalyst File Gen",    self.mfg,                 4, 5,
             "Creates CSV files for MetaboAnalyst"),
             
			("RNA-seq Deglist",      self.RNASeqDeg,                1, 6,
             "DESeq2 differential gene list generated from count matrix"),
            
            ("Gene -> Pathways",          self.gene2path,           2, 6,
             "Look up KEGG pathways associated with genes"),

            ("External Metrics",          self.externalOpt,         3, 6,
             "Compares clustering solutions with Rand-index, Adj. Rand-index, normalized mutual info, adj. mutual info"),

            ("VIP (MS/MS Comp.)",         self.progen,              4, 6,
             "Variable importance projections (VIP) MS/MS comparison"),
             
            ("Help/Documentation",        self.helpOut,             0, 7,
             "Opens the help wiki"),
        ]

        # Group by row and center each row inside the 4 button columns (1..4).
        by_row = defaultdict(list)
        for label, cmd, col, row, tip in buttons:
            by_row[row].append((label, cmd, tip))

        button_cols = [1, 2, 3, 4]
        for row in sorted(by_row.keys()):
            row_buttons = by_row[row]
            n = len(row_buttons)
            # Center rows within the 4 equal-width button columns (1..4).
            # For n=3 this yields columns 1..3 (left-biased), which matches the
            # visual centering when column 0 remains a smaller left padding column.
            start_idx = (len(button_cols) - n) // 2
            for i, (label, cmd, tip) in enumerate(row_buttons):
                b = ttk.Button(self, text=label, style="RW.TButton", command=cmd)
                b.grid(column=button_cols[start_idx + i], row=row, sticky=(N, S, E, W), padx=5, pady=5)
                CreateToolTip(b, tip)

    #  home  (now just delegates to _show_home)
    def home(self):
        self._show_home()

    #  sub-page scaffold
    def _sub_page(self, title):
        """Clear frame, apply sub-page grid, draw title + home button."""
        self._clear()
        self._configure_sub_grid()

        spacer = self._SUB_SPACER_ROW
        footer = self._SUB_FOOTER_ROW
        # Push footer to the bottom of the window when it is tall enough.
        self.rowconfigure(spacer, weight=1)

        ttk.Label(self, text=title,
                  font=("TkHeadingFont", 28, "bold")).grid(
            column=0, row=0, columnspan=5, pady=15, sticky=(N,))

        ttk.Button(self, text="Return to Home",
                   command=self.home).grid(
            column=0, row=footer, columnspan=5,
            pady=(12, 20), padx=60, sticky=(S, E, W))

    #  preproces
    def preprocess(self):
        filename     = filedialog.askopenfilename()
        if not filename:
            return
        metab_data   = GB.fileCheck(file=filename)
        metab_data_c = GB.fileCheck(file=filename)
        metab_data_c = metab_data.drop(0, axis=0)
        columns      = list(metab_data_c.columns)
        metab_data_c = metab_data_c.drop(columns[0], axis=1)
        metab_data_c = metab_data_c.drop(columns[-1], axis=1)
        metab_data_c = metab_data_c.to_numpy()
        data, toDelete = GB.dataCheck(metab_data_c)
        for i in range(len(toDelete)):
            toDelete[i] += 1
        metab_data = metab_data.drop(toDelete)
        metab_data.to_excel("pre_processed_data.xlsx", index=False)

    #  createClustergram  
    def createClustergram(self):
        self._sub_page("Clustergram Input")
    
        def linkageOutput(*args):
            self.sel_link = linkLB.curselection()
            if not self.sel_link: return
            curLink = linkageList[self.sel_link[0]]
            if curLink == 'ward':
                _refill(distLB, [distList[0]])
            else:
                _refill(distLB, distList)
    
        def submit(*args):
            try:
                norm_sel      = normList[normLB.curselection()[0]]
                scale_sel     = scaleList[scaleLB.curselection()[0]]
                transform_sel = transformList[transformLB.curselection()[0]]
                link_sel      = linkageList[linkLB.curselection()[0]]
                dist_sel      = distList[distLB.curselection()[0]]
                color_sel     = colorList[colorLB.curselection()[0]]
                groupOrd      = [g for g in grpVar.get().split(' ') if g.strip()]
                if not groupOrd:
                    groupOrd = ['1']
    
                if norm_sel == 'Normalize':
                    GU.createClustergram(1,
                                         linkFunc  = link_sel,
                                         distMet   = dist_sel,
                                         cmap      = color_sel,
                                         colOrder  = groupOrd,
                                         transform = transform_sel,
                                         scale     = 'NormStand')
                else:
                    n = 2 if len(groupOrd) > 1 else 0
                    GU.createClustergram(n,
                                         linkFunc  = link_sel,
                                         distMet   = dist_sel,
                                         cmap      = color_sel,
                                         colOrder  = groupOrd,
                                         transform = transform_sel,
                                         scale     = scale_sel)
            except Exception as e:
                import traceback
                traceback.print_exc()
    
        normList      = config.normList
        scaleList     = config.scaleList
        transformList = config.transformList
        linkageList   = config.linkageList
        distList      = config.distList
        colorList     = config.colorList
    
        for col, lbl in [(1,"Linkage"),(2,"Distance"),(3,"Color-Map"),
                         (4,"Transform"),(5,"Scale"),(6,"Normalize")]:
            ttk.Label(self, text=lbl,
                      font=("TkHeadingFont", 11)).grid(column=col, row=1, pady=3)
    
        linkLB      = _lb(self, col=1)
        distLB      = _lb(self, col=2)
        colorLB     = _lb(self, col=3)
        transformLB = _lb(self, col=4)
        scaleLB     = _lb(self, col=5)
        normLB      = _lb(self, col=6)
    
        _refill(linkLB,      linkageList)
        _refill(distLB,      distList)
        _refill(colorLB,     colorList)
        _refill(transformLB, transformList)
        _refill(scaleLB,     scaleList)
        _refill(normLB,      normList)
    
        self.update_idletasks()
    
        linkLB.bind('<Double-1>', linkageOutput)

        self._wiki_buttons_cmap_trans_scale_under_menus(3)

        ttk.Label(self, text="Group order (space separated):",
                  font=("TkHeadingFont", 11)).grid(
            column=1, row=4, columnspan=3, sticky='w', padx=10)
        grpVar = tk.StringVar()
        ttk.Entry(self, textvariable=grpVar).grid(
            column=1, row=5, columnspan=3, sticky='ew', padx=10)
        ttk.Button(self, text="Submit", command=submit).grid(
            column=4, row=6, padx=10, sticky='ew')

    #  clusterSelection                                                    
    def clusterSelection(self):
        self._sub_page("Cluster Selection")

        def linkageOutput(*args):
            global selection
            selection = distLB.curselection()
            curLink   = linkageList[selection[0]]
            _refill(sampleLB, [distList[0]] if curLink == 'ward' else distList)

        def colorMap(*args):
            global selection1
            selection1 = sampleLB.curselection()
            _refill(colorLB, colorList)

        def dataTransform(*args):
            global selection2
            selection2 = colorLB.curselection()
            _refill(transformLB, transformList)

        def dataScale(*args):
            global selection3
            selection3 = transformLB.curselection()
            _refill(scaleLB, scaleList)

        def dataNorm(*args):
            global selection4
            selection4 = scaleLB.curselection()
            _refill(normLB, normList)
            submit_btn.grid(column=4, row=7, columnspan=3,
                            padx=10, sticky=(E, W))

        def submit(*args):
            norm      = normList[normLB.curselection()[0]]
            dist      = distList[selection1[0]]
            link      = linkageList[selection[0]]
            color     = colorList[selection2[0]]
            transform = transformList[selection3[0]]
            scale     = scaleList[selection4[0]]
            config.colorNum = 0
            groupOrd  = grpVar.get().split(',')
            if norm == 'Normalize':
                GU.selectClusters(link, dist, 1, colOrder=groupOrd,
                                  transform=transform, scale='NormStand',
                                  cmap=color)
            else:
                GU.selectClusters(link, dist,
                                  transform=transform, scale=scale, cmap=color)

        linkageList   = config.linkageList
        distList      = config.distList
        colorList     = config.colorList
        transformList = config.transformList
        scaleList     = config.scaleList
        normList      = config.normList

        for col, lbl in [(1,"Linkage"),(2,"Distance"),(3,"Color-Map"),
                         (4,"Transform"),(5,"Scale"),(6,"Normalize?")]:
            ttk.Label(self, text=lbl,
                      font=("TkHeadingFont", 11)).grid(column=col, row=1, pady=3)

        distLB = _lb(self, 2, col=1)
        sampleLB = _lb(self, 2, col=2)
        colorLB = _lb(self, 2, col=3)
        transformLB = _lb(self, 2, col=4)
        scaleLB = _lb(self, 2, col=5)
        normLB = _lb(self, 2, col=6)

        _refill(distLB, linkageList)

        distLB.bind('<Double-1>',      linkageOutput)
        sampleLB.bind('<Double-1>',    colorMap)
        colorLB.bind('<Double-1>',     dataTransform)
        transformLB.bind('<Double-1>', dataScale)
        scaleLB.bind('<Double-1>',     dataNorm)

        self._wiki_buttons_cmap_trans_scale_under_menus(3)

        ttk.Label(self, text="Group order (normalization col first, separate by a space)",
                  font=("TkHeadingFont", 11)).grid(
            column=1, row=4, columnspan=3, sticky=W, padx=10)
        grpVar = tk.StringVar()
        ttk.Entry(self, textvariable=grpVar).grid(
            column=1, row=5, columnspan=3, sticky=(E, W), padx=10)
        submit_btn = ttk.Button(self, text="Submit", command=submit)

    #  medians
    def medians(self):
        self._sub_page("Group Medians")
        var1 = IntVar()
        ttk.Checkbutton(self, text="Remove Zeros?", variable=var1).grid(
            column=1, row=2, columnspan=3, pady=10)
        ttk.Button(self, text="Select file",
                   command=lambda: GU.groupMedians(rmZeros=var1.get())).grid(
            column=1, row=3, columnspan=3, sticky=(E, W), padx=60)

    #  linkages
    def linkages(self):
        self._sub_page("Linkage Comparison")

        def distFunc(*args):
            global distanceMet
            distanceMet = self.dist.get()
            values = [1,2,3,4] if distanceMet == 'euclidean' else [1,2,3]
            num_comps = StringVar()
            self.numComps = ttk.Combobox(self, values=values, textvariable=num_comps)
            self.numComps.bind('<<ComboboxSelected>>', numCompsFunc)
            self.numComps.grid(column=2, row=2, padx=5, pady=5)

        def numCompsFunc(*args):
            global numberComps
            numberComps = self.numComps.get()
            linkOpts = ['ward-single','ward-complete','ward-average',
                        'single-complete','single-average','complete-average',
                        'ward-single-complete','ward-single-average',
                        'ward-complete-average','single-complete-average',
                        'ward-single-complete-average']
            if numberComps == '1':
                value = (['ward','single','complete','average']
                         if distanceMet == 'euclidean'
                         else ['single','complete','average'])
            elif numberComps == '2':
                value = linkOpts[0:6] if distanceMet == 'euclidean' else linkOpts[3:6]
            elif numberComps == '3':
                value = linkOpts[6:10] if distanceMet == 'euclidean' else [linkOpts[9]]
            else:
                value = [linkOpts[10]]
            lv = StringVar()
            self.linkage = ttk.Combobox(self, values=value, textvariable=lv)
            self.linkage.bind('<<ComboboxSelected>>', linkageComp)
            self.linkage.grid(column=3, row=2, padx=5, pady=5)

        def linkageComp(*args):
            global linkList
            linkList  = []
            sel       = self.linkage.get()
            first     = 0
            for i, c in enumerate(sel):
                if c == '-':
                    linkList.append(sel[first:i]); first = i+1
                elif i == len(sel)-1:
                    linkList.append(sel[first:])
            if not linkList:
                linkList = [sel]
            _refill(transformLB, transformList)

        def dataScale(*args):
            global dataTrans
            dataTrans = transformList[transformLB.curselection()[0]]
            _refill(scaleLB, scaleList)

        def submit(*args):
            scale_ = scaleList[scaleLB.curselection()[0]]
            file   = filedialog.askopenfilename()
            GU.linkageComparison(file, numberComps, linkList,
                                 distanceMet, dataTrans, scale_)

        transformList = config.transformList
        scaleList     = config.scaleList
        distList      = config.distList

        ttk.Label(self, text="Distance measure",
                  font=("TkHeadingFont",12)).grid(column=1, row=1, pady=5)
        ttk.Label(self, text="# comparisons",
                  font=("TkHeadingFont",12)).grid(column=2, row=1, pady=5)
        ttk.Label(self, text="Linkage functions",
                  font=("TkHeadingFont",12)).grid(column=3, row=1, pady=5)
        ttk.Label(self, text="Transform",
                  font=("TkHeadingFont",12)).grid(column=1, row=3, pady=5)
        ttk.Label(self, text="Scale",
                  font=("TkHeadingFont",12)).grid(column=2, row=3, pady=5)

        distances = StringVar()
        self.dist = ttk.Combobox(self, values=distList, textvariable=distances)
        self.dist.bind('<<ComboboxSelected>>', distFunc)
        self.dist.grid(column=1, row=2, padx=5, pady=5)

        transformLB = _lb(self, 4, col=1)
        scaleLB     = _lb(self, 4, col=2)
        transformLB.bind('<Double-1>', dataScale)

        self._wiki_buttons_transform_scale(5, start_col=1)
        ttk.Button(self, text="Submit", command=submit).grid(
            column=1, row=6, columnspan=3, sticky=(E,W), padx=60, pady=5)

    #  heatmapAnalyses

    def heatmapAnalyses(self):
        self._sub_page("Heatmap Input")

        def linkageOutput(*args):
            global selection
            selection = distLB.curselection()
            _refill(sampleLB,
                    [distList[0]] if linkageList[selection[0]]=='ward' else distList)

        def colorMap(*args):
            global selection1
            selection1 = sampleLB.curselection()
            _refill(colorLB, colorList)

        def dataTransform(*args):
            global selection2
            selection2 = colorLB.curselection()
            _refill(transformLB, transformList)

        def dataScale(*args):
            global selection3
            selection3 = transformLB.curselection()
            _refill(scaleLB, scaleList)

        def dataNorm(*args):
            global selection4
            selection4 = scaleLB.curselection()
            _refill(normLB, normList)
            submit_btn.grid(column=4, row=7, columnspan=3,
                            padx=10, sticky=(E,W))

        def submit(*args):
            norm      = normList[normLB.curselection()[0]]
            dist      = distList[selection1[0]]
            link      = linkageList[selection[0]]
            color     = colorList[selection2[0]]
            transform = transformList[selection3[0]]
            scale     = scaleList[selection4[0]]
            groupOrd  = grpVar.get().split(',')
            if norm == 'Normalize':
                GU.heatmapAnalysis(link, dist, color, 1,
                                   colOrder=groupOrd,
                                   transform=transform, scale='NormStand')
            else:
                GU.heatmapAnalysis(link, dist, color, 0,
                                   transform=transform, scale=scale)

        linkageList   = config.linkageList
        distList      = config.distList
        colorList     = config.colorList
        transformList = config.transformList
        scaleList     = config.scaleList
        normList      = config.normList

        for col, lbl in [(1,"Linkage"),(2,"Distance"),(3,"Color-Map"),
                         (4,"Transform"),(5,"Scale"),(6,"Normalize?")]:
            ttk.Label(self, text=lbl,
                      font=("TkHeadingFont",11)).grid(column=col, row=1, pady=3)

        # place each Listbox directly under its header label to keep
        # the inputs visually aligned with their titles
        distLB      = _lb(self, 2, col=1)  # Linkage
        sampleLB    = _lb(self, 2, col=2)  # Distance
        colorLB     = _lb(self, 2, col=3)  # Color-Map
        transformLB = _lb(self, 2, col=4)  # Transform
        scaleLB     = _lb(self, 2, col=5)  # Scale
        normLB      = _lb(self, 2, col=6)  # Normalize?

        _refill(distLB, linkageList)

        distLB.bind('<Double-1>',      linkageOutput)
        sampleLB.bind('<Double-1>',    colorMap)
        colorLB.bind('<Double-1>',     dataTransform)
        transformLB.bind('<Double-1>', dataScale)
        scaleLB.bind('<Double-1>',     dataNorm)

        self._wiki_buttons_cmap_trans_scale_under_menus(3)

        ttk.Label(self, text="Group order (normalization col first, separate by a space)",
                  font=("TkHeadingFont",11)).grid(
            column=1, row=4, columnspan=3, sticky=W, padx=10)
        grpVar = tk.StringVar()
        ttk.Entry(self, textvariable=grpVar).grid(
            column=1, row=5, columnspan=3, sticky=(E,W), padx=10)
        submit_btn = ttk.Button(self, text="Submit", command=submit)

    #  thin-delegate methods                                               
    def compound(self):   GU.compoundMatchUp(typeFile='all')
    def integrity(self):
        f = filedialog.askopenfilename()
        if f: GU.dataIntegrity(f)
    def mstF(self):
        numClust = self._run_with_busy(
            "Working... running MST optimization.",
            lambda: GU.MST(self, func='ensemble')
        )
        if numClust is not None:
            GU.ensembleClustering(optNum=int(numClust))


    def P2P(self):        GU.peaksToPathways()
    def helpOut(self):    webbrowser.open('https://github.com/hisl6802/ClusteringToolbox/wiki')
    def progen(self):     GU.progenesis()
    def gene2path(self):  GU.geneToPathway()
    def RNASeqDeg(self): GU.RNASeqDeg()
    
    #  ECCO
    def ensemble(self):
        self._sub_page("Ensemble Clustering with Optimization")

        def ensembleCMap(*args):
            global cMapSel
            cMapSel = colorList[colorMapLB.curselection()[0]]
            _refill(transformLB, transformList)

        def ensembleTransform(*args):
            global transSel
            transSel = transformList[transformLB.curselection()[0]]
            _refill(scaleLB, scaleList)

        def ensembleScale(*args):
            global scaleSel
            scaleSel = scaleList[scaleLB.curselection()[0]]
            submit_btn.grid(column=1, row=5, columnspan=3,
                            sticky=(E,W), padx=60, pady=5)

        colorList     = config.colorList
        transformList = config.transformList
        scaleList     = config.scaleList

        for col, lbl in [(1,"ColorMap"),(2,"Transform"),(3,"Scale")]:
            ttk.Label(self, text=lbl,
                      font=("TkHeadingFont",16)).grid(column=col, row=1, pady=5)

        colorMapLB  = _lb(self, 2, col=1)
        transformLB = _lb(self, 2, col=2)
        scaleLB     = _lb(self, 2, col=3)
        _refill(colorMapLB, colorList)

        colorMapLB.bind('<Double-1>',  ensembleCMap)
        transformLB.bind('<Double-1>', ensembleTransform)
        scaleLB.bind('<Double-1>',     ensembleScale)

        self._wiki_buttons_color_transform_scale(3)

        # checkbox to run Peaks->Pathways and Compound MatchUp cascade
        self.run_p2p_cascade = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            self,
            text="Also generate PeaksToPathways files and prompt for Mummichog results",
            variable=self.run_p2p_cascade
        ).grid(column=1, row=4, columnspan=3, pady=5, padx=5, sticky=(E, W))

        def run_ensemble_and_cascade():
            sel = colorMapLB.curselection()
            col_map = colorList[sel[0]] if sel else 'viridis'
            result = GU.ensembleClusteringFullOpt(
                transform=transSel, scale=scaleSel, colMap=col_map,
            )
            try:
                opt_k, out_dir = result
            except Exception:
                # older behavior: ignore if return value not unpackable
                return

            if not self.run_p2p_cascade.get():
                return

            # attempt to run PeaksToPathways using the same raw file
            raw_file = getattr(config, 'last_ensemble_file', None)
            if raw_file and out_dir:
                GU.peaksToPathways(raw_file=raw_file, clusters_dir=out_dir)

            # prompt once for a folder of Mummichog results and run Compound MatchUp
            messagebox.showinfo(
                title="Mummichog results",
                message=(
                    "Select the folder containing your Mummichog matched-compound CSV files\n"
                    "(for example, the output folder from MummiBot)."
                ),
            )
            mummi_dir = filedialog.askdirectory()
            if mummi_dir:
                GU.compoundMatchUp(typeFile='all', file=mummi_dir)

        submit_btn = ttk.Button(
            self, text="Submit",
            command=run_ensemble_and_cascade)

    def cooccHeatmap(self):
        self._sub_page('CoOcc Heatmap')

        coocc_var = tk.StringVar(value='')
        data_var = tk.StringVar(value='')
        k_var = tk.StringVar(value='4')

        transformList = config.transformList
        scaleList = config.scaleList
        co_ts = {'transform': 'None', 'scale': 'None'}

        ttk.Label(
            self,
            text='Number of clusters (1–10):',
            font=('TkHeadingFont', 12),
        ).grid(column=1, row=1, sticky=W, padx=10, pady=6)

        ttk.Combobox(
            self,
            textvariable=k_var,
            values=[str(i) for i in range(1, 11)],
            width=6,
            state='readonly',
        ).grid(column=2, row=1, sticky=W, padx=5, pady=6)

        cmap_list = list(dict.fromkeys(config.colorList))
        if 'coolwarm' in cmap_list:
            cmap_list.remove('coolwarm')
        cmap_list = ['coolwarm'] + cmap_list
        cmap_var = tk.StringVar(value='coolwarm')

        ttk.Label(
            self,
            text='Colormap:',
            font=('TkHeadingFont', 12),
        ).grid(column=1, row=2, sticky=W, padx=10, pady=4)
        ttk.Combobox(
            self,
            textvariable=cmap_var,
            values=cmap_list,
            width=16,
            state='readonly',
        ).grid(column=2, row=2, sticky=W, padx=5, pady=4)
        ttk.Button(
            self,
            text='ColorMap Options',
            command=lambda: webbrowser.open(self._URL_COLORMAP),
        ).grid(column=3, row=2, padx=5, pady=4, sticky=W)

        for col, lbl in [(1, 'Transform'), (2, 'Scale')]:
            ttk.Label(
                self, text=lbl, font=('TkHeadingFont', 12),
            ).grid(column=col, row=3, pady=(6, 2), padx=10, sticky=W)

        def on_coocc_trans(*_):
            if not trans_lb.curselection():
                return
            co_ts['transform'] = transformList[trans_lb.curselection()[0]]
            _refill(scale_lb, scaleList)

        def on_coocc_scale(*_):
            if not scale_lb.curselection():
                return
            co_ts['scale'] = scaleList[scale_lb.curselection()[0]]

        trans_lb = _lb(self, 4, col=1, height=6)
        scale_lb = _lb(self, 4, col=2, height=6)
        _refill(trans_lb, transformList)
        trans_lb.bind('<Double-1>', on_coocc_trans)
        scale_lb.bind('<Double-1>', on_coocc_scale)

        self._wiki_buttons_transform_scale(5, start_col=1)

        cluster_cols_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            self,
            text='Cluster columns by group (column dendrogram)',
            variable=cluster_cols_var,
        ).grid(column=1, row=6, columnspan=3, sticky=W, padx=10, pady=2)

        ttk.Label(
            self,
            text='EnsembleCoOcc.xlsx:',
            font=('TkHeadingFont', 12),
        ).grid(column=1, row=7, sticky=W, padx=10, pady=4)

        def browse_coocc():
            p = filedialog.askopenfilename(
                title='Select EnsembleCoOcc.xlsx',
                filetypes=[
                    ('Excel', '*.xlsx'),
                    ('All files', '*.*'),
                ],
            )
            if p:
                coocc_var.set(p)

        ttk.Entry(
            self, textvariable=coocc_var, width=56,
        ).grid(column=1, row=8, columnspan=2, sticky=(E, W), padx=10, pady=2)
        ttk.Button(self, text='Browse…', command=browse_coocc).grid(
            column=3, row=8, padx=5, pady=2, sticky=W,
        )

        ttk.Label(
            self,
            text='Group medians matrix (.xlsx, row 1 = group names):',
            font=('TkHeadingFont', 12),
        ).grid(column=1, row=9, sticky=W, padx=10, pady=4)

        def browse_data():
            p = filedialog.askopenfilename(
                title='Select group medians matrix (e.g. *_matrix_medians.xlsx)',
                filetypes=[
                    ('Excel', '*.xlsx'),
                    ('CSV', '*.csv'),
                    ('All files', '*.*'),
                ],
            )
            if p:
                data_var.set(p)

        ttk.Entry(
            self, textvariable=data_var, width=56,
        ).grid(column=1, row=10, columnspan=2, sticky=(E, W), padx=10, pady=2)
        ttk.Button(self, text='Browse…', command=browse_data).grid(
            column=3, row=10, padx=5, pady=2, sticky=W,
        )

        def run_coocc():
            self._run_with_busy(
                'Working… CoOcc heatmaps.',
                lambda: GU.cooccHeatmap(
                    k_var.get(),
                    coocc_var.get(),
                    data_var.get(),
                    cmap=cmap_var.get(),
                    col_cluster=cluster_cols_var.get(),
                    transform=co_ts['transform'],
                    scale=co_ts['scale'],
                ),
                status_row=11,
            )

        ttk.Button(self, text='Run', command=run_coocc).grid(
            column=1, row=12, columnspan=2, sticky=(E, W), padx=60, pady=15,
        )

    #  MST
    def mst(self):
        self._sub_page("Validation Clusters (5 metrics)")

        def transType(*args):
            global transform
            transform = transformList[transLB.curselection()[0]]
            _refill(scaleLB, scaleList)

        def scaleType(*args):
            global scale
            scale = scaleList[scaleLB.curselection()[0]]
            _refill(valLB, valList)

        def valType(*args):
            sel = valList[valLB.curselection()[0]]
            self._run_with_busy(
                f"Working... running {sel} validation.",
                lambda: GU.MST(self, transform=transform, scale=scale, func=sel)
            )

        transformList = config.transformList
        scaleList     = config.scaleList
        valList       = ('k-means based','DBI','Dunn','PBM','Silhouette','All')

        for col, lbl in [(1,"Transform"),(2,"Scale"),(3,"Validation")]:
            ttk.Label(self, text=lbl,
                      font=("TkHeadingFont",14)).grid(column=col, row=1, pady=5)

        transLB = _lb(self, 2, col=1)
        scaleLB = _lb(self, 2, col=2)
        valLB   = _lb(self, 2, col=3)
        _refill(transLB, transformList)

        transLB.bind('<Double-1>', transType)
        scaleLB.bind('<Double-1>', scaleType)
        valLB.bind('<Double-1>',   valType)

        self._wiki_buttons_transform_scale(3)

    #  genSelClustFig                                                      
    def genSelClustFig(self):
        self._sub_page("Selected Clusters Figure")
        colorList = config.colorList
        ttk.Label(self, text="Select ColorMap:",
                  font=("TkHeadingFont",14)).grid(column=1, row=1, columnspan=3)
        cmap_lb = _lb(self, 2, col=1, colspan=3)
        _refill(cmap_lb, colorList)
        cmap_lb.bind('<Double-1>',
                     lambda *_: GB.createHeatmapFig(
                         clMap=colorList[cmap_lb.curselection()[0]]))
        ttk.Button(
            self,
            text='ColorMap Options',
            command=lambda: webbrowser.open(self._URL_COLORMAP),
        ).grid(column=2, row=3, padx=5, pady=5, sticky=(E, W))


    #  anovaHeatMap
    def anovaHeatMap(self):
        self._sub_page("ANOVA Heatmap")

        def transformAH(*args):
            global cMap
            cMap = colorList[colorMapLB.curselection()[0]]
            _refill(transformLB, transformList)

        def scaleAHMP(*args):
            global transform
            transform = transformList[transformLB.curselection()[0]]
            _refill(scaleLB, scaleList)
            submit_btn.grid(column=1, row=6, columnspan=3,
                            sticky=(E,W), padx=60, pady=5)

        colorList     = config.colorList
        transformList = config.transformList
        scaleList     = config.scaleList

        for col, lbl in [(1,"ColorMap"),(2,"Transform"),(3,"Scale")]:
            ttk.Label(self, text=lbl,
                      font=("TkHeadingFont",16)).grid(column=col, row=1, pady=5)

        colorMapLB  = _lb(self, 2, col=1)
        transformLB = _lb(self, 2, col=2)
        scaleLB     = _lb(self, 2, col=3)
        _refill(colorMapLB, colorList)

        colorMapLB.bind('<Double-1>',  transformAH)
        transformLB.bind('<Double-1>', scaleAHMP)

        self._wiki_buttons_color_transform_scale(3)

        submit_btn = ttk.Button(
            self, text="Submit",
            command=lambda: GU.anovaHM(
                transform=transform,
                scale=scaleList[scaleLB.curselection()[0]],
                cMap=cMap))

    #  enzymeLookUp                                                     
    def enzymeLookUp(self):
        self._sub_page("Enzyme Look Up")
        numHM = tuple(range(1, 15))
        ttk.Label(self, text="Number of Heatmap sheets:",
                  font=("TkHeadingFont",14)).grid(column=1, row=1, columnspan=3)
        lb = _lb(self, 2, col=1, colspan=3)
        _refill(lb, numHM)
        lb.bind('<Double-1>',
                lambda *_: GU.enzymeLookUp(
                    numSheets=numHM[lb.curselection()[0]]))

    #  CIsTtest
    def CIsTtest(self):
        self._sub_page("Confidence Intervals from t-tests")

        def submitNumSamps(*args):
            global numSampsSel
            numSampsSel = numSamps[numSampsList.curselection()[0]]
            submit_btn.grid(column=1, row=5, columnspan=3,
                            sticky=(E,W), padx=60, pady=5)

        numSamps = tuple(range(1, 101))
        ttk.Label(self, text="Samples per group:",
                  font=("TkHeadingFont",13)).grid(column=1, row=1, columnspan=3)
        numSampsList = _lb(self, 2, col=1, colspan=3)
        _refill(numSampsList, numSamps)
        numSampsList.bind('<Double-1>', submitNumSamps)

        ttk.Label(self, text="Confidence level (e.g. 95):",
                  font=("TkHeadingFont",13)).grid(column=1, row=3, columnspan=3)
        self.conf = tk.StringVar()
        ttk.Entry(self, textvariable=self.conf).grid(
            column=1, row=4, columnspan=3, sticky=(E,W), padx=60)

        submit_btn = ttk.Button(
            self, text="Submit",
            command=lambda: GU.confidenceIntervals(
                numSampsSel, confidenceLevel=self.conf.get()))

    #  normality check                                                         
    def normalityC(self):
        self._sub_page("Normality Check")

        def transSelected(*args):
            global curTrans
            curTrans = transformList[transformLB.curselection()[0]]
            _refill(scaleLB, scaleList)
            submit_btn.grid(column=1, row=5, columnspan=3,
                            sticky=(E,W), padx=60, pady=5)

        transformList = config.transformList
        scaleList     = config.scaleList

        for col, lbl in [(1,"Transform"),(2,"Scale")]:
            ttk.Label(self, text=lbl,
                      font=("TkHeadingFont",14)).grid(column=col, row=1, pady=5)

        transformLB = _lb(self, 2, col=1)
        scaleLB     = _lb(self, 2, col=2)
        _refill(transformLB, transformList)
        transformLB.bind('<Double-1>', transSelected)

        self._wiki_buttons_transform_scale(3)

        submit_btn = ttk.Button(
            self, text="Submit",
            command=lambda: GU.normalityCheck(
                transform=curTrans,
                scale=scaleList[scaleLB.curselection()[0]]))

    #  mfg                                                                 
    def mfg(self):
        self._sub_page("Metaboanalyst File Generation")
        options = ["Multi", "Uni", "All"]
        ttk.Label(self, text="Select variant:",
                  font=("TkHeadingFont",14)).grid(column=1, row=1, columnspan=3)
        lb = _lb(self, 2, col=1, colspan=3)
        _refill(lb, options)

        def submitMFG(*args):
            selected = options[lb.curselection()[0]]
            file     = filedialog.askopenfilename()
            if file:
                GU.mfgUtil(fileName=file, variant=selected)

        lb.bind('<Double-1>', submitMFG)

    #  externalOpt                                                         
    def externalOpt(self):
        self._sub_page("External Criteria Evaluation")

        def transformation(*args):
            global transSelection
            transSelection = config.transformList[transformLB.curselection()[0]]
            _refill(scaleLB, config.scaleList)

        def scaling(*args):
            global scaleSelection
            scaleSelection = config.scaleList[scaleLB.curselection()[0]]
            _refill(metricsLB, config.metrics)

        def submitExternal(*args):
            sel = config.metrics[metricsLB.curselection()[0]]
            GU.externalCriteria(comp=sel,
                                trans=transSelection,
                                scale=scaleSelection)

        for col, lbl in [(1,"Transformation"),(2,"Scale"),(3,"Criteria")]:
            ttk.Label(self, text=lbl,
                      font=("TkHeadingFont",14)).grid(column=col, row=1, pady=5)

        transformLB = _lb(self, 2, col=1)
        scaleLB     = _lb(self, 2, col=2)
        metricsLB   = _lb(self, 2, col=3)
        _refill(transformLB, config.transformList)

        transformLB.bind('<Double-1>', transformation)
        scaleLB.bind('<Double-1>',     scaling)
        metricsLB.bind('<Double-1>',   submitExternal)

    #  monoVal                                                             
    def monoVal(self):
        self._sub_page("Mono-clustering Optimization")

        def dataTransform(*args):
            global selection2
            selection2 = optLB.curselection()
            _refill(transformLB, transformList)

        def dataScale(*args):
            global selection3
            selection3 = transformLB.curselection()
            _refill(scaleLB, scaleList)
            
        def linkageOutput(*args):
            global selection
            selection = distLB.curselection()
            _refill(sampleLB,
                    [distList[0]] if linkageList[selection[0]]=='ward' else distList)

        def optMet(*args):
            global selection1
            selection1 = sampleLB.curselection()
            _refill(optLB, optList)

        def upNumClust(*args):
            global selection4
            selection4 = scaleLB.curselection()
            _refill(clustersLB, minMetab)

        def submit(*args):
            sel5      = clustersLB.curselection()
            dist      = distList[selection1[0]]
            link      = linkageList[selection[0]]
            optimizer = optList[selection2[0]]
            transform = transformList[selection3[0]]
            scale     = scaleList[selection4[0]]
            numClust  = minMetab[sel5[0]]
            GU.validateMono(self, numClust, transform, scale,
                            link, dist, optimizer)

        linkageList   = config.linkageList
        distList      = config.distList
        optList       = config.optList
        transformList = config.transformList
        scaleList     = config.scaleList
        minMetab      = tuple(range(3, 101))

        for col, lbl in [(1,"Linkage"),(2,"Distance"),(3,"Opt. Metric"),
                         (4,"Transform"),(5,"Scale"),(6,"Max clusters")]:
            ttk.Label(self, text=lbl,
                      font=("TkHeadingFont",11)).grid(column=col, row=1, pady=3)

        distLB      = _lb(self, 2, col=1); sampleLB    = _lb(self, 2, col=2)
        optLB       = _lb(self, 2, col=3); transformLB = _lb(self, 2, col=4)
        scaleLB     = _lb(self, 2, col=5); clustersLB  = _lb(self, 2, col=6)

        _refill(distLB, linkageList)

        distLB.bind('<Double-1>',      linkageOutput)
        sampleLB.bind('<Double-1>',    optMet)
        optLB.bind('<Double-1>',       dataTransform)
        transformLB.bind('<Double-1>', dataScale)
        scaleLB.bind('<Double-1>',     upNumClust)

        self._wiki_buttons_transform_scale(3, start_col=4)

        ttk.Button(self, text="Submit", command=submit).grid(
            column=1, row=4, columnspan=6, sticky=(E,W), padx=60, pady=8)


#  module-level helpers                                                
def _refill(lb, items):
    lb.delete(0, tk.END)
    for i, item in enumerate(items):
        lb.insert(i, item)

def _lb(parent, *args, col=None, row=None, colspan=1, rowspan=1, height=8):
    """
    Create a Listbox with consistent defaults.

    Backwards compatible with the older helper signature _lb(parent, col, row=2)
    AND newer call-sites that pass (row, col=...) or (row, col=..., colspan=...).
    """
    if len(args) > 2:
        raise TypeError("_lb() takes at most 3 positional arguments (parent, col|row, row)")

    # Defaults
    _row = 2 if row is None else row
    _col = 1 if col is None else col

    if len(args) == 1:
        a0 = args[0]
        if col is not None and row is None:
            # _lb(parent, row, col=...)
            _row = a0
            _col = col
        elif col is None and row is None:
            # _lb(parent, col)
            _col = a0
            _row = 2
        else:
            # If row explicitly provided, treat positional as col
            _col = a0
    elif len(args) == 2:
        # _lb(parent, col, row)
        _col, _row = args

    lb = tk.Listbox(parent, height=height, exportselection=False)
    lb.grid(column=_col, row=_row, columnspan=colspan, rowspan=rowspan,
            sticky='nsew', padx=5, pady=5)
    return lb

#  entry point                                                         
if __name__ == '__main__':
    multiprocessing.freeze_support()
    root = tk.Tk()
    root.minsize(700, 500)
    root.title("Clustering Toolbox")
    app = JuneLabClusteringGUI(master=root)
    root.mainloop()