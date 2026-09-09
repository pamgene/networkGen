test_that("generate_networks_batch dispatches each task, merges extra_args, and carries meta through untouched", {
  future::plan(future::sequential)

  mock_generate <- function(x, ppi_network, write = TRUE) {
    list(x = x, ppi_network = ppi_network, write = write)
  }

  tasks <- list(
    list(args = list(x = 1), meta = list(condition = "A", role = "observed")),
    list(args = list(x = 2, write = FALSE), meta = list(condition = "A", role = "permutation", perm_index = 1)),
    list(args = list(x = 3), meta = list(condition = "B", role = "observed"))
  )

  out <- generate_networks_batch(
    tasks, generate_fn = mock_generate, ppi_network = "fake_ppi",
    extra_args = list(write = TRUE)
  )

  expect_length(out, 3)
  expect_equal(out[[1]]$result$x, 1)
  expect_equal(out[[1]]$result$write, TRUE) # from extra_args
  expect_equal(out[[2]]$result$write, FALSE) # task-level arg wins over extra_args
  expect_equal(out[[1]]$result$ppi_network, "fake_ppi") # shared across every task
  expect_equal(out[[2]]$meta, list(condition = "A", role = "permutation", perm_index = 1))
  expect_equal(out[[3]]$meta$condition, "B")
})

test_that("generate_networks_batch preserves task order and lets a task's result be NULL", {
  future::plan(future::sequential)

  mock_generate <- function(fail, ppi_network) if (fail) NULL else "ok"

  tasks <- list(
    list(args = list(fail = FALSE), meta = list(id = 1)),
    list(args = list(fail = TRUE), meta = list(id = 2)),
    list(args = list(fail = FALSE), meta = list(id = 3))
  )

  out <- generate_networks_batch(tasks, generate_fn = mock_generate, ppi_network = NULL)

  expect_equal(purrr::map_int(out, ~ .x$meta$id), c(1, 2, 3))
  expect_equal(out[[1]]$result, "ok")
  expect_null(out[[2]]$result)
  expect_equal(out[[3]]$result, "ok")
})
