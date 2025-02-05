AUTOMAKE_OPTIONS = foreign

#Put here all the workflows and associate python files here
WORKFLOWFILES = crires_classification.py crires_datasources.py crires_keywords.py crires_rules.py crires_science.py crires_standards.py crires_task_functions.py crires_wavelength.py crires_wkf.py __init__.py

#Put here all the recipe parameters associated with the workflow
PARAMETERFILES = crires_parameters.yaml

#This installs the workflow and the OCA rules in reflexaux directory
workflow_DATA = $(WORKFLOWFILES) ${PARAMETERFILES}

EXTRA_DIST = $(WORKFLOWFILES) ${PARAMETERFILES}

