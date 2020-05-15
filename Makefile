SLV= slv-v0.5/slv
CCI= cci-v0.0/cci

TARGETS= slv cci

.PHONY: all clean FORCE

all: $(TARGETS)

FORCE:

clean:
	cd $(dir $(SLV)) && $(MAKE) clean
	cd $(dir $(CCI)) && $(MAKE) clean
	rm -f $(TARGETS)

$(SLV): FORCE
	cd $(dir $@) && $(MAKE)

$(CCI): $(SLV) FORCE
	cd $(dir $@) && $(MAKE)

slv: $(SLV)
	cp $< $@

cci: $(CCI)
	cp $< $@
