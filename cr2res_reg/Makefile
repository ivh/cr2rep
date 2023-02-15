TESTS ?= $(wildcard *.yaml)
GEN ?= "eso-testgen.py"
GOLD ?= "gold"
RUNSCRIPTS := $(patsubst %.yaml,%.mk, $(TESTS))

export ESOREX_PLUGIN_DIR=$(CURDIR)/../recipes/.libs
all: validate $(RUNSCRIPTS)

validate: $(TESTS)
	if python -c 'import jsonschema' 2>/dev/null; then $(GEN) --validate $^; fi

#generate NRI/p_regtest compatible recipe configurations
create_sof: $(TESTS) validate
	$(GEN) --flatten  $(TESTS)

%.mk: %.yaml
	$(GEN) --forgive -o "new/$^" -g "$(GOLD)/$^" -r "$@" -t makefile "$^"

# run test scripts and generate or verify against a gold dataset
check: validate $(RUNSCRIPTS)
	@if test ! -d "${SOF_DATA}" ; then \
                echo "ERROR: The SOF_DATA environment variable was not set to point to the regression data." 1>&2 ; \
                exit 1 ; \
        fi

	set -e && for r in $(RUNSCRIPTS); do \
	    $(MAKE) -f $$r; \
	done

valgrind: validate
	set -e && for r in $(TESTS); do \
	    $(GEN) --valgrind -o valgrind/$$r -r $${r}-valgrind.mk -t makefile $$r; \
	    $(MAKE) -f $${r}-valgrind.mk; \
	done

# compare already existing test results with fitsdiff-folder
compare:
	set -e && for r in $(TESTS); do \
	  [ ! -e new/$$r ] && echo "new/$$r missing run check first" && continue; \
	  fitsdiff-folder.sh --late-fail $(GOLD)/$$r new/$$r \
	    --ignore-keywords='DATE,ESO PRO REC*,DATAMD5,CHECKSUM,DATASUM'; \
	done

# create pipeline regression report in xml and html formats
# the following environment variables must be set
# DATA_REF : with reference results of pipeline products created by make check
# PIPETOOLS_DIR : directory containing pipetools project
report:
	# test required env variables are set, else exit with error message
	@if test ! -d "${DATA_REF}" ; then \
                echo "ERROR: The DATA_REF environment variable must be set to point to the reference regression results data." 1>&2 ; \
                exit 1 ; \
        fi

	@if test ! -d "${PIPETOOLS_DIR}" ; then \
                echo "ERROR: The PIPETOOLS_DIR environment variable must be set to point to the directory where the SVN pipetools project is checkout." 1>&2 ; \
                exit 1 ; \
        fi

	# set additional required env variables
	export DATA_TEST=${PWD}/gold ; \
	export DATA_YAML_DIR=${PWD} ; \
	cd ${PIPETOOLS_DIR}/pipetools/production/pipeline_regtest_py/ ; \
	python -m pytest pipeline_reg_tests*.py -v --tb=short --junit-xml=regression_test_results.xml --html=regression_html_results.html --self-contained-html ; \
	cp regression_test_results.xml regression_html_results.html ${DATA_YAML_DIR} ; \
	cd -



#clean all generated files
clean:
	rm -rf new
	rm -f *.mk
	rm -f valgrind.sh
	rm -rf valgrind
	rm -rf gold_rrrecipe*.html
	rm -rf gold_rrrecipe_calib*.html
	rm -rf regression_html_results.html  regression_test_results.xml
	find sof -mindepth 1 -type d  | xargs rm -rf

# clean up all generated files and the reference dataset
distclean: clean
	rm -rf $(GOLD)

.PHONY: all check clean distclean nri validate
.NOTPARALLEL:
