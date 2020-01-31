
AUTOMAKE_OPTIONS = foreign

#Put here all the workflows. They must end with extension .in. They should
#be like a normal workflow but using @prefix@ for paths 
#to the OCA rules (grep for OCA in the xml file) or Python scripts
WORKFLOWS = cr2res.xml
#Put here the OCA rules used by the workflow
OCAWKF =  cr2res_wkf.oca
#Put Python scripts used by the workflows
PYTHONWKF = interactive_window.py

#This installs the workflow and the OCA rules in reflexaux directory
wkfextra_DATA = $(WORKFLOWS) $(OCAWKF) $(PYTHONWKF)

#Pack in the distribution the original workflow
EXTRA_DIST = $(WORKFLOWS).in $(OCAWKF) $(PYTHONWKF)

#This also installs the workflow in the common reflex directory
wkfcopy_DATA = $(WORKFLOWS)

