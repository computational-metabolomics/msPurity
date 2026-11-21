test_that("MSP comparisons preserve content across platform line endings", {
    lines <- c("NAME: example", "Num Peaks: 2", "100 10", "200 20", "")
    expected <- "NAME: example\nNum Peaks: 2\n100 10\n200 20\n"
    for (ending in c(LF = "\n", CRLF = "\r\n", CRCRLF = "\r\r\n")) {
        f <- tempfile(fileext = ".msp")
        writeBin(charToRaw(paste(lines, collapse = ending)), f)
        expect_identical(.msp_text(f), expected,
                         info = paste("line ending", encodeString(ending)))
    }
})
