'use strict';

/* Controllers */

var viperAppControllers = angular.module('viperAppControllers', []);

viperAppControllers.controller('FullReportCtrl', ['$scope', '$http',
  function($scope, $http) {
    $http.get('app/model/report.json').success(function(data) {
      $scope.report = data;
    });
  }]);

viperAppControllers.controller('SectionReportCtrl', ['$scope', '$http', '$routeParams',
  function($scope, $http, $routeParams) {
    $http.get('app/model/' + $routeParams.sectionName + '.json').success(function(data) {
      $scope.report = data;
    });
  }]);
