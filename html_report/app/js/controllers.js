'use strict';

/* Controllers */

var rnaAppControllers = angular.module('rnaAppControllers', []);

rnaAppControllers.controller('FullReportCtrl', ['$scope', '$http',
  function($scope, $http) {
    $http.get('model/report.json').success(function(data) {
      $scope.report = data;
    });

    $scope.orderProp = 'sectionId';
  }]);

rnaAppControllers.controller('SectionReportCtrl', ['$scope', '$http',
  function($scope, $http) {
    $http.get('model/' + $routeParams.sectionId + '.json').success(function(data) {
      $scope.report = data;
    });
  }]);
