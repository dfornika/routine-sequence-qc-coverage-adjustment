(ns coverage-adjust.core
  (:require [clojure.java.io :as io]
            [clojure.tools.cli :refer [parse-opts]]
            [clojure.data.csv :as csv])
  (gen-class))

(def cli-options
  [["-b" "--basic-qc BASIC_QC"]
   ["-h" "--Help"]])

(defn -main [& args]
  (def opts (parse-opts args cli-options))
  #_(def opts {:options {}
               :arguments []
               :summary "  -c, --config CONFIG\n  -h, --Help"
               :errors nil})
  (pprint opts)
  )

(comment
  )