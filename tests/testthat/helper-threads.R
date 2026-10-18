# Package checks run small fixtures with one native thread. Individual tests
# may temporarily change it and restore this value.
getFromNamespace("set_scNPM_threads", "scNPM")(1L)
