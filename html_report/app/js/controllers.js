'use strict';

/* Controllers */

var viperAppControllers = angular.module('viperAppControllers', []);

viperAppControllers.controller('FullReportCtrl', ['$scope', '$http',
  function($scope, $http) {
    $http.get('app/model/report.json').success(function(data) {
      $scope.report = data;
    });
  }]);

viperAppControllers.controller('SectionReportCtrl', ['$scope', '$http', '$location',
  function($scope, $http, $location) {
    var cur_path = $location.path();
    $http.get('/app/model' + cur_path + '.json').success(function(data) {
      $scope.report = data;
    });
  }]);
