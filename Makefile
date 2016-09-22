# This file is part of DBSCAN.
# https://github.com/josl/DBSCAN

# Licensed under the MIT license:
# http://www.opensource.org/licenses/MIT-license
# Copyright (c) 2016, Jose L. Bellod Cisneros <bellod.cisneros@gmail.com>

# lists all available targets
list:
	@sh -c "$(MAKE) -p no_targets__ | awk -F':' '/^[a-zA-Z0-9][^\$$#\/\\t=]*:([^=]|$$)/ {split(\$$1,A,/ /);for(i in A)print A[i]}' | grep -v '__\$$' | grep -v 'make\[1\]' | grep -v 'Makefile' | sort"
# required for list
no_targets__:

# install all dependencies (do not forget to create a virtualenv first)
setup:
	@pip install -U -e .\[tests\]

# test your application (tests in the tests/ directory)
test: mongo_test unit

unit:
	@coverage run --branch `which nosetests` -vv --with-yanc -s tests/
	@coverage report -m --fail-under=80

# show coverage in html format
coverage-html: unit
	@coverage html

# get a mongodb instance up (localhost:3333)
mongo: kill_mongo
	@mkdir -p /tmp/DBSCAN/mongodata && mongod --dbpath /tmp/DBSCAN/mongodata --logpath /tmp/DBSCAN/mongolog --port 3333 --quiet &

# kill this mongodb instance (localhost:3333)
kill_mongo:
	@-ps aux | egrep -i 'mongod.+3333' | egrep -v egrep | awk '{ print $$2 }' | xargs kill -9

# clear all data in this mongodb instance (localhost: 3333)
clear_mongo:
	@rm -rf /tmp/DBSCAN && mkdir -p /tmp/DBSCAN/mongodata

# get a mongodb instance up for your unit tests (localhost:3334)
mongo_test: kill_mongo_test
	@rm -rf /tmp/DBSCAN/mongotestdata && mkdir -p /tmp/DBSCAN/mongotestdata
	@mongod --dbpath /tmp/DBSCAN/mongotestdata --logpath /tmp/DBSCAN/mongotestlog --port 3334 --quiet --fork
	@echo 'waiting for mongo...'
	@until mongo --port 3334 --eval "quit()"; do sleep 0.25; done > /dev/null 2> /dev/null

# kill the test mongodb instance (localhost: 3334)
kill_mongo_test:
	@-ps aux | egrep -i 'mongod.+3334' | egrep -v egrep | awk '{ print $$2 }' | xargs kill -9

# run tests against all supported python versions
tox:
	@tox

#docs:
	#@cd DBSCAN/docs && make html && open _build/html/index.html
