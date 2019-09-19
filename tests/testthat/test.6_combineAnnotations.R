context ("checking combine annotations based functions")

test_that("checking combine annotations based functions", {
  print ("\n")
  print("########################################################")
  print("## Checking combine annotations based functions       ##")
  print("########################################################")
  metfrag_resultPth <- system.file("extdata", "tests", "external_annotations", "metfrag.tsv", package="msPurity")
  sirius_csi_resultPth <- system.file("extdata","tests", "external_annotations", "sirus_csifingerid.tsv", package="msPurity")
  probmetab_resultPth <- system.file("extdata","tests", "external_annotations", "probmetab.tsv", package="msPurity")
  ms1_lookup_resultPth <- system.file("extdata","tests", "external_annotations", "beams.tsv", package="msPurity")
  sm_resultPth <- system.file("extdata","tests", "sm", "spectralMatching_result.sqlite", package="msPurity")
  compoundDbPth <- system.file("extdata","tests", "db", "metab_compound_subset.sqlite", package="msPurity")


  #localCompoundDbPth <- system.file("extdata", "tests", "db", "compounds_18July2019_0319.sqlite", package="msPurity")
  weights=list('sm'=0.3,
               'metfrag'=0.2,
               'sirius_csifingerid'=0.2,
               'probmetab'=0, # ignore probmetab results now that the package is no longer maintained
               'ms1_lookup'=0.05,
               'biosim'=0.25
  )
  td <- tempdir()
  sm_resultPthCopy <- file.path(td, 'sm_result_tmp41.sqlite')
  file.copy(sm_resultPth, sm_resultPthCopy)

  sqlitePth <- sm_resultPthCopy

  combined <- combineAnnotations(sqlitePth,
                                 compoundDbPth=compoundDbPth,
                                 metfrag_resultPth = metfrag_resultPth,
                                 sirius_csi_resultPth = sirius_csi_resultPth,
                                 probmetab_resultPth = probmetab_resultPth,
                                 ms1_lookup_resultPth = ms1_lookup_resultPth,
                                 weights = weights,
                                 ms1_lookup_dbSource='hmdb',
                                 ms1_lookup_checkAdducts=FALSE,
                                 ms1_lookup_keepAdducts=c('[M+H]+', '[M-H]-')
                                 )
  # likely to change to only providing limited testing
  expect_equal(combined$rank,   c(1,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,
                                  5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,7,7,7,7,8,8,8,9,9,9,9,9,9,9,9,9,9,
                                  9,9,9,9,9,9,9,9,10,11,11,12,12,12,12,13,13,13,13,13,13,13,13,14,
                                  14,15,16,16,16,16,17,18,18,19,19,20,21,1,1,1,2,3,1))

  # sql command to get reduced list of compounds from database
  # DELETE FROM pubchem WHERE inchikey NOT IN (
  # 'AGPKZVBTJJNPAG-UHFFFAOYSA-N',
  # 'ARJPPNFIEQKVBB-UHFFFAOYSA-N',
  # 'DHMQDGOQFOQNFH-UHFFFAOYSA-N',
  # 'FTWHFXMUJQRNBK-UHFFFAOYSA-N',
  # 'HNJBEVLQSNELDL-UHFFFAOYSA-N',
  # 'IAZDPXIOMUYVGZ-UHFFFAOYSA-N',
  # 'LQNUZADURLCDLV-UHFFFAOYSA-N',
  # 'MWFMGBPGAXYFAR-UHFFFAOYSA-N',
  # 'OAICVXFJPJFONN-UHFFFAOYSA-N',
  # 'ONIBWKKTOPOVIA-AZXPZELESA-N',
  # 'ONIBWKKTOPOVIA-BBKVAIMGSA-N',
  # 'ONIBWKKTOPOVIA-BFEYZEMLSA-N',
  # 'ONIBWKKTOPOVIA-BQNZLWMBSA-N',
  # 'ONIBWKKTOPOVIA-BWEKNVLDSA-N',
  # 'ONIBWKKTOPOVIA-BYPYZUCNSA-M',
  # 'ONIBWKKTOPOVIA-BYPYZUCNSA-N',
  # 'ONIBWKKTOPOVIA-BYPYZUCNSA-O',
  # 'ONIBWKKTOPOVIA-DGLZGNQDSA-N',
  # 'ONIBWKKTOPOVIA-DHWWSWJHSA-N',
  # 'ONIBWKKTOPOVIA-DYCDLGHISA-N',
  # 'ONIBWKKTOPOVIA-FMMSUUDPSA-N',
  # 'ONIBWKKTOPOVIA-FQXRHUMFSA-N',
  # 'ONIBWKKTOPOVIA-FVSUZMELSA-N',
  # 'ONIBWKKTOPOVIA-GDYNRSRVSA-N',
  # 'ONIBWKKTOPOVIA-GIZBTRSZSA-N',
  # 'ONIBWKKTOPOVIA-GTTLGWSSSA-N',
  # 'ONIBWKKTOPOVIA-HOSYLAQJSA-N',
  # 'ONIBWKKTOPOVIA-HRCSTIOUSA-N',
  # 'ONIBWKKTOPOVIA-IDEBNGHGSA-N',
  # 'ONIBWKKTOPOVIA-IJDGHPMYSA-N',
  # 'ONIBWKKTOPOVIA-ITEPJMEFSA-N',
  # 'ONIBWKKTOPOVIA-IXBOUXNVSA-N',
  # 'ONIBWKKTOPOVIA-JGTYJTGKSA-N',
  # 'ONIBWKKTOPOVIA-JQZHSJCGSA-N',
  # 'ONIBWKKTOPOVIA-JRGPAWSWSA-N',
  # 'ONIBWKKTOPOVIA-JUYFNQAYSA-N',
  # 'ONIBWKKTOPOVIA-KICNZHNUSA-N',
  # 'ONIBWKKTOPOVIA-KIZNEYSQSA-N',
  # 'ONIBWKKTOPOVIA-KRGXXEMVSA-N',
  # 'ONIBWKKTOPOVIA-LCBCNYHZSA-N',
  # 'ONIBWKKTOPOVIA-LTDLRDEHSA-N',
  # 'ONIBWKKTOPOVIA-OFHXNXSMSA-N',
  # 'ONIBWKKTOPOVIA-OHMILCFJSA-N',
  # 'ONIBWKKTOPOVIA-OHXCBXKRSA-N',
  # 'ONIBWKKTOPOVIA-OUBTZVSYSA-N',
  # 'ONIBWKKTOPOVIA-OXZJVQSUSA-N',
  # 'ONIBWKKTOPOVIA-OYYOGNGZSA-N',
  # 'ONIBWKKTOPOVIA-OZJLVGIFSA-N',
  # 'ONIBWKKTOPOVIA-PEPZRWDSSA-N',
  # 'ONIBWKKTOPOVIA-PEPZRWDSSA-O',
  # 'ONIBWKKTOPOVIA-PJDIGCCDSA-N',
  # 'ONIBWKKTOPOVIA-PMCDIQQISA-N',
  # 'ONIBWKKTOPOVIA-PPHVBSBQSA-N',
  # 'ONIBWKKTOPOVIA-PTQBSOBMSA-N',
  # 'ONIBWKKTOPOVIA-QDNHWIQGSA-N',
  # 'ONIBWKKTOPOVIA-QOOYSLLASA-N',
  # 'ONIBWKKTOPOVIA-QQTTXZKSSA-N',
  # 'ONIBWKKTOPOVIA-QRTGCQPVSA-N',
  # 'ONIBWKKTOPOVIA-QSBWZAQZSA-N',
  # 'ONIBWKKTOPOVIA-QYKNYGDISA-N',
  # 'ONIBWKKTOPOVIA-RHRFEJLCSA-N',
  # 'ONIBWKKTOPOVIA-RVQWGROCSA-N',
  # 'ONIBWKKTOPOVIA-SCSAIBSYSA-M',
  # 'ONIBWKKTOPOVIA-SCSAIBSYSA-N',
  # 'ONIBWKKTOPOVIA-SCSAIBSYSA-O',
  # 'ONIBWKKTOPOVIA-SETFFETJSA-N',
  # 'ONIBWKKTOPOVIA-SFIIULIVSA-N',
  # 'ONIBWKKTOPOVIA-SLJODTCNSA-N',
  # 'ONIBWKKTOPOVIA-SRVZCWNMSA-N',
  # 'ONIBWKKTOPOVIA-TXZHAAMZSA-N',
  # 'ONIBWKKTOPOVIA-UBBKNGMPSA-N',
  # 'ONIBWKKTOPOVIA-UHFFFAOYSA-M',
  # 'ONIBWKKTOPOVIA-UHFFFAOYSA-N',
  # 'ONIBWKKTOPOVIA-UHFFFAOYSA-O',
  # 'ONIBWKKTOPOVIA-VMGGCIAMSA-N',
  # 'ONIBWKKTOPOVIA-VNEWRNQKSA-N',
  # 'ONIBWKKTOPOVIA-VNFUZYCESA-N',
  # 'ONIBWKKTOPOVIA-VSDVZINISA-N',
  # 'ONIBWKKTOPOVIA-WAPLMMNXSA-N',
  # 'ONIBWKKTOPOVIA-WGGUOBTBSA-N',
  # 'ONIBWKKTOPOVIA-WRHOHBQNSA-N',
  # 'ONIBWKKTOPOVIA-XAFSXMPTSA-N',
  # 'ONIBWKKTOPOVIA-XQRIOMIJSA-N',
  # 'ONIBWKKTOPOVIA-YIKKUXJFSA-N',
  # 'ONIBWKKTOPOVIA-YZRVCBOHSA-N',
  # 'ONIBWKKTOPOVIA-ZYXKZDFDSA-N',
  # 'PVNIIMVLHYAWGP-UHFFFAOYSA-N',
  # 'QQZWEECEMNQSTG-UHFFFAOYSA-N',
  # 'RRUDCFGSUDOHDG-UHFFFAOYSA-N',
  # 'RWRDLPDLKQPQOW-CNRUNOGKSA-N',
  # 'RWRDLPDLKQPQOW-DICFDUPASA-N',
  # 'RWRDLPDLKQPQOW-DYCDLGHISA-N',
  # 'RWRDLPDLKQPQOW-HJOWPTDZSA-N',
  # 'RWRDLPDLKQPQOW-KHORGVISSA-N',
  # 'RWRDLPDLKQPQOW-KLRAWXKOSA-N',
  # 'RWRDLPDLKQPQOW-LBPDFUHNSA-N',
  # 'RWRDLPDLKQPQOW-LNLMKGTHSA-N',
  # 'RWRDLPDLKQPQOW-MICDWDOJSA-N',
  # 'RWRDLPDLKQPQOW-MNYXATJNSA-N',
  # 'RWRDLPDLKQPQOW-QAOQSSEZSA-N',
  # 'RWRDLPDLKQPQOW-SMZGMGDZSA-N',
  # 'RWRDLPDLKQPQOW-SVYQBANQSA-N',
  # 'RWRDLPDLKQPQOW-UHFFFAOYSA-N',
  # 'RWRDLPDLKQPQOW-UHFFFAOYSA-O',
  # 'RWRDLPDLKQPQOW-WFVSFCRTSA-N',
  # 'RWRDLPDLKQPQOW-YZRHJBSPSA-N',
  # 'SIOXPEMLGUPBBT-UHFFFAOYSA-N',
  # 'TWBYWOBDOCUKOW-UHFFFAOYSA-N',
  # 'VULIHENHKGDFAB-UHFFFAOYSA-N',
  # 'XYLMUPLGERFSHI-UHFFFAOYSA-N')


})
