	.PHONEY: all save clean

	all:
	for d in $$(ls -d ./Practical_Programming/*); do \
		$(MAKE) -C $$d; \
	done;
	for d in $$(ls -d ./Numerical_Methods/*); do \
		$(MAKE) -C $$d; \
	done;
	$(MAKE) -C ./Exam;


	save :
		git add -A
		git commit -m "Auto save";
		git push origin master;


	clean :
	for d in $$(ls -d ./Practical_Programming/*); do \
		$(MAKE) -C $$d clean; \
	done;
	for d in $$(ls -d ./Numerical_Methods/*); do \
		$(MAKE) -C $$d clean; \
	done;
	$(MAKE) -C ./Exam clean;

